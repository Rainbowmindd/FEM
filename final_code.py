import math
import numpy as np


class Node():
    def __init__(self,x,y,BC=False):
        self.x=x
        self.y=y
        self.BC=BC #true or false , 1 or 0
        
class Element():
    def __init__(self,ID,Jakobian=None):
        
        self.ID=ID
        self.Jakobian=Jakobian
        self.H = [[0.0 for _ in range(4)] for _ in range(4)]  
        self.Hbc = [[0.0 for _ in range(4)] for _ in range(4)]
        self.P_local = [0.0 for _ in range(4)]
        self.P=[0.0 for _ in range(4)] #P for each element
        self.C_local=[[0.0 for _ in range(4)] for _ in range(4)]
        self.C=[[0.0 for _ in range(4)] for _ in range(4)]
        
class SOE():
    def __init__(self,Nn):
        self.H_global=[[0.0 for _ in range(Nn)] for _ in range(Nn)]
        self.C_global= [[0.0 for _ in range(Nn)] for _ in range(Nn)]
        self.P_global=[0.0 for _ in range(Nn)]
       
class GaussIntegration:
    def __init__(self, npc):
        if npc == 4:  
            self.pc = [-0.5773502691896257, 0.5773502691896257]
            self.w = [1.0, 1.0]
        elif npc == 9:
            self.pc = [-0.7745966692414834, 0.0, 0.7745966692414834]
            self.w = [5/9, 8/9, 5/9]
        elif npc == 16: 
            self.pc = [-0.8611363115940526, -0.3399810435848563, 0.3399810435848563, 0.8611363115940526]
            self.w = [0.3478548451374538, 0.6521451548625461, 0.6521451548625461, 0.3478548451374538]
        else:
            raise ValueError("Error")    
    
class Grid():
    def __init__(self,nN,nE):
        self.nN=nN #number of nodes
        self.nE=nE #number of elements
        self.nodes=[]
        self.elements=[]
        
    def display_nodes(self):
        print("\nWspolrzedne wezlow i BC: \n")
        for i, node in enumerate(self.nodes):
            print(f"Wezel {i+1}: ({node.x}, {node.y}, BC: {node.BC})")
            
    def display_elements(self):
        print("\nId wezlow poszczegolnych elementow: \n")
        for i, element in enumerate(self.elements):
            print(f"Element {i+1}: Wezly {element.ID}")
 
class ElementUniv: 
    def __init__(self,npc):
        self.dN_dksi=[] 
        self.dN_deta=[] 
        self.N_wart_funkcji=[]
        
        gaussIntegration = GaussIntegration(npc)
        points = gaussIntegration.pc 
        
        for ksi in points:
            for eta in points:
                self.dN_dksi.append([-0.25 * (1 - eta),
                                     0.25 * (1 - eta),
                                     0.25 * (1 + eta),
                                     -0.25 * (1 + eta)])

                self.dN_deta.append([-0.25 * (1 - ksi),
                                     -0.25 * (1 + ksi),
                                     0.25 * (1 + ksi),
                                     0.25 * (1 - ksi)])
        
                N1 = 0.25 * (1 - ksi) * (1 - eta)
                N2 = 0.25 * (1 + ksi) * (1 - eta)
                N3 = 0.25 * (1 + ksi) * (1 + eta)
                N4 = 0.25 * (1 - ksi) * (1 + eta)
                
                self.N_wart_funkcji.append([N1, N2, N3, N4])

        #self.printTabs()
        
    def printTabs(self):
        print("dN/dksi table: ")
        for i in range(len(self.dN_dksi)):
            print(self.dN_dksi[i])
        print("\ndN/deta table:  ")
        for i in range(len(self.dN_deta)):
            print(self.dN_deta[i])  
        print("\nWartosci funkcji ksztaltu w punktach calkowania:")
        for i in range(len(self.N_wart_funkcji)):
            print(self.N_wart_funkcji[i])         
class Jakobian:
    def __init__(self, nodeElement, elementUniv, npc):
        self.J = [[0.0 for _ in range(4)] for _ in range(npc)]
        self.J1 = [[0.0 for _ in range(4)] for _ in range(npc)]
        self.detJ = []
        for l in range(npc):
            for i in range(4):
                self.J[l][0] += elementUniv.dN_dksi[l][i] * nodeElement[i].x
                self.J[l][1] += elementUniv.dN_dksi[l][i] * nodeElement[i].y
                self.J[l][2] += elementUniv.dN_deta[l][i] * nodeElement[i].x
                self.J[l][3] += elementUniv.dN_deta[l][i] * nodeElement[i].y

            
            epsilon=0.00001
            # elementy !=0
            self.J[l] = [0 if abs(val) < epsilon else val for val in self.J[l]]

          
            detJ_val = self.J[l][0] * self.J[l][3] - (self.J[l][1] * self.J[l][2])
            self.detJ.append(detJ_val)

            if detJ_val > 0:
                self.J1[l][0] = (1 / detJ_val) * self.J[l][3]
                self.J1[l][1] = (1 / detJ_val) * -self.J[l][1]
                self.J1[l][2] = (1 / detJ_val) * -self.J[l][2]
                self.J1[l][3] = (1 / detJ_val) * self.J[l][0]
            else:
                self.J1[l] = [0, 0, 0, 0]

def agregation(self,element,local_H,C_local):
    global_H=[[0 for _ in range(4)] for _ in range(4)]
    global_C=[[0 for _ in range(4)] for _ in range(4)] 
    
    ids_of_node=element.ID
    for i in range(4):
        for j in range(4):
            global_i=ids_of_node[i]-1 
            global_j=ids_of_node[j]-1
            global_H[global_i][global_j] +=local_H[i][j]
    for i in range(4):
        for j in range(4):
            global_i=ids_of_node[i]-1 
            global_j=ids_of_node[j]-1
            global_C[global_i][global_j]+=C_local[i][j]



class GlobalData():
    def __init__(self, lines):
       #dane symulacji
            self.SimulationTime=int(lines[0].split()[1])
            self.SimulationStepTime=int(lines[1].split()[1])
            self.Conductivity=int(lines[2].split()[1])
            self.Alfa=int(lines[3].split()[1])
            self.Tot=int(lines[4].split()[1])
            self.InitialTemp=int(lines[5].split()[1])
            self.Density=int(lines[6].split()[1])
            self.SpecificHeat=int(lines[7].split()[1])
            self.nN=int(lines[8].split()[2]) 
            self.nE=int(lines[9].split()[2])
      
def load_data_from_file(lines,grid):
            #Szukanie *Node
            i=0 
            while i<len(lines) and lines[i].strip() !="*Node":
                i+=1
            i+=1  
            if i >= len(lines):
                raise ValueError("Brak wezlow w pliku")
            
            for j in range(grid.nN):
                temp=lines[i].strip().split(",")
                x,y=float(temp[1]),float(temp[2])
                grid.nodes.append(Node(x, y))
                i+=1
            #Szukanie *Element
            while i<len(lines) and lines[i].strip() !="*Element, type=DC2D4":
                i+=1
            i+=1
            if i >= len(lines):
                raise ValueError("Brak elementow w pliku")
            
         
            for j in range(grid.nE):
                temp=lines[i].strip().split(",")
                nodeID=list(map(int,temp[1:]))
                grid.elements.append(Element(nodeID))
                i+=1
               
               #Szukanie *BC
            while i < len(lines) and lines[i].strip() != "*BC":
                i += 1
            i += 1  
            if i >= len(lines):
                raise ValueError("Brak BC w pliku")

            temp = lines[i].strip().split(",")
            boundary_conditions = list(map(int, temp))
            for bc in boundary_conditions:
                grid.nodes[bc - 1].BC = True #ustawienie flagi dla odpowiednich wezlow
            
            while i < len(lines) and "InitialTemp" not in lines[i]:
                i += 1
    
            if i < len(lines) and "InitialTemp" in lines[i]:
                initial_temp = float(lines[i].split("=")[1].strip())
                for node in grid.nodes:
                    node.temp = initial_temp
                  
                     
def calcH(jakobian, elementUniv,element_nodes,surface, k,alfa, npc):
    dN_dx = [[0 for _ in range(4)] for _ in range(npc)]
    dN_dy = [[0 for _ in range(4)] for _ in range(npc)]
    for m, j1 in enumerate(jakobian.J1):
        for x in range(4):
            dN_dx[m][x] = j1[0] * elementUniv.dN_dksi[m][x] + j1[1] * elementUniv.dN_deta[m][x]
            dN_dy[m][x] = j1[2] * elementUniv.dN_dksi[m][x] + j1[3] * elementUniv.dN_deta[m][x]

    # Mnozenie transponownych i zwyklych
    HpcX = [[[0 for _ in range(4)] for _ in range(4)] for _ in range(npc)]
    HpcY = [[[0 for _ in range(4)] for _ in range(4)] for _ in range(npc)]
    
    for integrPoint in range(npc):
        for x in range(4):
            for y in range(4):
                HpcX[integrPoint][x][y] = dN_dx[integrPoint][x] * dN_dx[integrPoint][y]
                HpcY[integrPoint][x][y] = dN_dy[integrPoint][x] * dN_dy[integrPoint][y]


    Hpc = [[[0 for _ in range(4)] for _ in range(4)] for _ in range(npc)]
    
    for integrPoint in range(npc):
        for y in range(4):
            for x in range(4):
                Hpc[integrPoint][y][x] = k * (HpcX[integrPoint][y][x] + HpcY[integrPoint][y][x]) * jakobian.detJ[integrPoint]
    gaussIntegration = GaussIntegration(npc)
    weights = gaussIntegration.w
    n = int(math.sqrt(npc))
        
        
    H = [[0 for _ in range(4)] for _ in range(4)]
    for i in range(npc):
        w1 = i % n 
        w2 = i//n 
        for x in range(4):
            for y in range(4):
                H[x][y] += Hpc[i][x][y] * weights[w1] * weights[w2]

    
    Hbc = calcHbc(surface, element_nodes, npc, alfa)

    for x in range(4):
        for y in range(4):
            H[x][y] += Hbc[x][y]
    return H


class Surface:
    def __init__(self, npc):
        self.N = []
        gaussIntegration = GaussIntegration(npc)
        points=gaussIntegration.pc
        weights=gaussIntegration.w
        
        for ksi in points:  #dla 1d wzor to N=(1-ksi)/2
            self.N.append([
                [(1 - ksi) / 2, (1 + ksi) / 2, 0, 0],  #gorny obrocony-> dolny bok
                [0, (1 - ksi) / 2, (1 + ksi) / 2, 0],  #odbity lewy-> prawy bok
                [0, 0, (1 + ksi) / 2, (1 - ksi) / 2],  #dolny obrocony -> gorny bok
                [(1 - ksi) / 2, 0, 0, (1 + ksi) / 2]   #odbity prawy-> lewy bok
            ])
               
def calcHbc(surface, element_nodes, npc, alfa):
    gaussIntegration = GaussIntegration(npc)
    Hbc = [[0.0 for _ in range(4)] for _ in range(4)]
    n = int(math.sqrt(npc))
    
    for side in range(4):
        if element_nodes[side].BC and element_nodes[(side + 1) % 4].BC:
            # Oblicz długość boku
            x1, y1 = element_nodes[side].x, element_nodes[side].y
            x2, y2 = element_nodes[(side + 1) % 4].x, element_nodes[(side + 1) % 4].y
            length = math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)
            detJ = length / 2.0
            
            for point in range(n):
                N = surface.N[point][side]
                weight = gaussIntegration.w[point]
            
                for i in range(4):
                    for j in range(4):
                        Hbc[i][j] += alfa * N[i] * N[j] * weight * detJ
    
    return Hbc

def calcP_Local(surface, element_nodes, npc, alfa):
    gaussIntegration = GaussIntegration(npc)
    
    P_local = [0.0 for _ in range(4)]
    n = int(math.sqrt(npc))
    
    for side in range(4):  # for every side
        if element_nodes[side].BC and element_nodes[(side + 1) % 4].BC:  # Check if the edge has boundary condition nodes
            x1, y1 = element_nodes[side].x, element_nodes[side].y
            x2, y2 = element_nodes[(side + 1) % 4].x, element_nodes[(side + 1) % 4].y
            length = math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)
            detJ = length / 2.0
            
            for point in range(n):
                N = surface.N[point][side]
                weight = gaussIntegration.w[point]   
                
                for i in range(4):
                    P_local[i] += alfa * N[i] * gaussIntegration.w[point] * detJ * 1200

    return P_local


def calc_global_P(grid, surface, npc, alfa):
    global_P = [0.0 for _ in range(len(grid.nodes))]  

   
    for element in grid.elements:
        element_nodes = [grid.nodes[id-1] for id in element.ID]  
        
        P_local = [0.0 for _ in range(4)]
        P_local = calcP_Local(surface, element_nodes, npc, alfa)
       
        for i in range(4):
            global_P[element.ID[i] - 1] += P_local[i]  

    return global_P

def calcC(element_nodes, jakobian, npc, c, ro, elementUniv):
    C_matrix = [[0.0 for _ in range(4)] for _ in range(4)]
    
    gaussIntegration = GaussIntegration(npc)
    weights = gaussIntegration.w
    n = int(math.sqrt(npc)) 
    
    for i in range(npc):
        w1 = i % n  # index dla wag ksi
        w2 = i // n  # index dla wag eta
        
        for x in range(4):
            for y in range(4):
                C_matrix[x][y] += (elementUniv.N_wart_funkcji[i][x] * 
                                 elementUniv.N_wart_funkcji[i][y] * 
                                 c * ro * jakobian.detJ[i] *
                                 weights[w1] * weights[w2])
    
    return C_matrix
       
def solve_temperature(global_H, global_P):
    """
   Uklad rownan H * t + P = 0
    """
    global_H = np.array(global_H)
    global_P = np.array(global_P)
    
    # Zmiana znaku w wektorze P
   # global_P = -global_P
    
    try:
        temperature = np.linalg.solve(global_H, global_P) #dekompozycja LU z pivotingiem
        print('\nRozwiazanie ukladu rownan:')
        print(temperature)
    except np.linalg.LinAlgError as e:
        print("\n")
        return None
    
    return temperature       

def time_solution(grid,data,surface,elem_univ,npc):
    
    soe=SOE(data.nN)
    
    #temperatures 
    t_0=np.full(data.nN,data.InitialTemp)
    t_1=np.zeros(data.nN)
    
    
    results=[]
    
    #time step loop
    t=0
    while t<data.SimulationTime:
        
        soe.H_global=[[0.0 for _ in range(data.nN)] for _ in range(data.nN)]
        soe.C_global=[[0.0 for _ in range(data.nN)] for _ in range(data.nN)]
        soe.P_global=[0.0 for _ in range(data.nN)]
      
        for element in grid.elements:
            element_nodes=[grid.nodes[id-1] for id in element.ID]
            
            jakobian=Jakobian(element_nodes,elem_univ,npc)
            H_local = calcH(jakobian, elem_univ, element_nodes, surface, 
                          data.Conductivity, data.Alfa, npc)
            C_local = calcC(element_nodes, jakobian, npc, 
                          data.SpecificHeat, data.Density, elem_univ)
            P_local = calcP_Local(surface, element_nodes, npc, data.Alfa)
            
            
            #aggregation
            for i in range(4):
                for j in range(4):
                    global_i = element.ID[i] - 1
                    global_j = element.ID[j] - 1
                    soe.H_global[global_i][global_j] += H_local[i][j]
                    soe.C_global[global_i][global_j] += C_local[i][j]
                soe.P_global[element.ID[i] - 1] += P_local[i]
        
        #list->array
        H=np.array(soe.H_global)
        C=np.array(soe.C_global)
        P=np.array(soe.P_global)
        
        #final: [H + C/dt]t1 = P + C/dt * t0
        dt=data.SimulationStepTime
        C_dt=C/dt
        A=H+C_dt
        b=P+np.dot(C_dt,t_0)
        
     
        try:
            t_1=np.linalg.solve(A,b)
            print(f"\nTemperatures at time {t + dt}:")
            print(t_1)
            results.append((np.min(t_1),np.max(t_1)))
        except np.linalg.LinAlgError as e:
            print("error")
            break
        
        t_0=t_1.copy()
        t+=dt
        
    print("\n Wyniki:")
    print("=" * 40)
    print(f"{'Time (s)':10} {'Min (°C)':20} {'Max (°C)':20}")
    print("-" * 40)

    for i, (min_temp, max_temp) in enumerate(results):
        time = (i + 1) * data.SimulationStepTime
        print(f"{time:<10.1f} {min_temp:<20.6f} {max_temp:<15.6f}")

    print("=" * 40)

    return t_1

def run_simulation(filename):
   
    with open(filename, 'r') as file:
        lines = file.readlines()
    
    data = GlobalData(lines)
    grid = Grid(data.nN, data.nE)
    load_data_from_file(lines, grid)
    
    #How many npc: 
    # npc=4
    npc=16
    # npc = 9 
    elem_univ = ElementUniv(npc)
    surface = Surface(npc)
    
 
    final_temperatures = time_solution(grid, data, surface, elem_univ, npc)
    
    return final_temperatures
          
 
#---Pliki----------------------------------
# final_temps = run_simulation('Test1_4_4.txt')
# final_temps = run_simulation('Test2_4_4_MixGrid.txt')
final_temps = run_simulation('Test3_31_31_kwadrat.txt')




  