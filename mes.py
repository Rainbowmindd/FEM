import math

class Node():
    def __init__(self,x,y,BC=False):
        self.x=x
        self.y=y
        self.BC=BC #true or false , 1 or 0
        
class Element():
    def __init__(self,ID,Jakobian=None): #domyslna wartosc dla jakobianu
        
        self.ID=ID
        #lab3:
        self.Jakobian=Jakobian
        #Hbc lab:
        self.H = [[0.0 for _ in range(4)] for _ in range(4)]  
        self.Hbc = [[0.0 for _ in range(4)] for _ in range(4)]
        
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
 
class ElementUniv: #element uniwersalny,niezalezny od siatki
    def __init__(self,npc):
        #dN_dksi and dN_deta (pochodne funkcji ksztaltu)
        self.dN_dksi=[] #npc x 4(nodes)
        self.dN_deta=[] 
        
        #Integration points
        points = [-0.7745966692414834, 0, 0.7745966692414834]  # Ksi, Eta 
       # weights = [5 / 9, 8 / 9, 5 / 9]

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

        self.printTabs()
        
    def printTabs(self):
        print("dN/dksi table: ")
        for i in range(len(self.dN_dksi)):
            print(self.dN_dksi[i])
        print("\ndN/deta table:  ")
        for i in range(len(self.dN_deta)):
            print(self.dN_deta[i])           
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
            # Upewnienie się, że elementy są różne od zera
            self.J[l] = [0 if abs(val) < epsilon else val for val in self.J[l]]

            # Obliczanie wyznacznika
            detJ_val = self.J[l][0] * self.J[l][3] - (self.J[l][1] * self.J[l][2])
            self.detJ.append(detJ_val)

            # Obliczanie odwrotności Jacobiana
            if detJ_val > 0:
                self.J1[l][0] = (1 / detJ_val) * self.J[l][3]
                self.J1[l][1] = (1 / detJ_val) * -self.J[l][1]
                self.J1[l][2] = (1 / detJ_val) * -self.J[l][2]
                self.J1[l][3] = (1 / detJ_val) * self.J[l][0]
            else:
                self.J1[l] = [0, 0, 0, 0]

class GlobalH:
    def __init__(self,size):
            self.globalH = [[0 for _ in range(size)] for _ in range(size)]   


def agregation(self,element,local_H):
    ids_of_node=element.ID
    for i in range(4):
        for j in range(4):
            global_i=ids_of_node[i]-1 
            global_j=ids_of_node[j]-1
            global_H[global_i][global_j] +=local_H[i][j]
            #global_H[global_i][global_j] = global_H[global_i][global_j] + local_H[i][j] +Hbc[i][j]



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
            
            #liczba wezlow i elementow
            self.nN=int(lines[8].split()[2]) 
            self.nE=int(lines[9].split()[2])
        
    
        
def load_data_from_file(lines,grid):
            #Szukanie *Node
            i=0 
            while i<len(lines) and lines[i].strip() !="*Node":
                i+=1
            i+=1  #przechodzenie do linijki z danymi wezlow
            if i >= len(lines):
                raise ValueError("Brak wezlow w pliku")
            #read nodes
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
            
            #read elements
            for j in range(grid.nE):
                temp=lines[i].strip().split(",")
                nodeID=list(map(int,temp[1:]))
                grid.elements.append(Element(nodeID))
                i+=1
               
               #Szukanie *BC
            while i < len(lines) and lines[i].strip() != "*BC":
                i += 1
            i += 1  # Przechodzenie do linii z danymi BC
            if i >= len(lines):
                raise ValueError("Brak BC w pliku")


            temp = lines[i].strip().split(",")
            boundary_conditions = list(map(int, temp))
            for bc in boundary_conditions:
                grid.nodes[bc - 1].BC = True #ustawienie flagi dla odpowiednich wezlow
                  
                     
#LAB 4 matrix H----------------------------
def calcH(jakobian, elementUniv, k, npc):
    dN_dx = [[0 for _ in range(4)] for _ in range(npc)]
    dN_dy = [[0 for _ in range(4)] for _ in range(npc)]
    for m, j1 in enumerate(jakobian.J1):
        for x in range(4):
            dN_dx[m][x] = j1[0] * elementUniv.dN_dksi[m][x] + j1[1] * elementUniv.dN_deta[m][x]
            dN_dy[m][x] = j1[2] * elementUniv.dN_dksi[m][x] + j1[3] * elementUniv.dN_deta[m][x]

    # print("\ndN/dx: ")
    # for i in range(len(dN_dx)):
    #     print(dN_dx[i])
    # print("\ndN/dy: ")
    # for i in range(len(dN_dy)):
    #     print(dN_dy[i])

    # Mnozenie transponownych i zwyklych
    HpcX = [[[0 for _ in range(4)] for _ in range(4)] for _ in range(npc)]
    HpcY = [[[0 for _ in range(4)] for _ in range(4)] for _ in range(npc)]
    
    for integrPoint in range(npc):
        for x in range(4):
            for y in range(4):
                HpcX[integrPoint][x][y] = dN_dx[integrPoint][x] * dN_dx[integrPoint][y]
                HpcY[integrPoint][x][y] = dN_dy[integrPoint][x] * dN_dy[integrPoint][y]

    # Obliczanie Hpc dla każdego punktu całkowania
    Hpc = [[[0 for _ in range(4)] for _ in range(4)] for _ in range(npc)]
    for integrPoint in range(npc):
        for y in range(4):
            for x in range(4):
                Hpc[integrPoint][y][x] = k * (HpcX[integrPoint][y][x] + HpcY[integrPoint][y][x]) * jakobian.detJ[integrPoint]

    # for i in range(len(Hpc)):
    #     print("Hpc", i + 1)
    #     for j in range(len(Hpc[i])):
    #         print(Hpc[i][j])
    #     print()

    # Obliczenie H
    weights = [5/9, 8/9, 5/9, 5/9, 8/9, 5/9, 5/9, 8/9, 5/9]  # Wagi dla 9 punktów całkowania
    H = [[0 for _ in range(4)] for _ in range(4)]
    for i in range(len(Hpc)):
        w1 = i % 3 #3x3 kwadratura gaussa przechodzimy przez indeksy 0 1 2 (ksi)
        w2 = (i // 3) % 3 #eta wynik=nr wiersza siatki
        for x in range(4):
            for y in range(4):
                H[x][y] += Hpc[i][x][y] * weights[w1] * weights[w2]

    print("Macierz H:")
    for row in H:
        print(row)
        
    return H

class Surface:
    def __init__(self,npc):
        self.N=[]
        gaussIntegration=GaussIntegration(npc)
        
        for ksi in gaussIntegration.pc:
            self.N.append([(1 - ksi) / 2, (1 + ksi) / 2, 0, 0])  #gorny i dolny bok
            self.N.append([0, (1 - ksi) / 2, (1 + ksi) / 2, 0])  #prawy bok
            self.N.append([0, 0, (1 + ksi) / 2, (1 - ksi) / 2])  #gorny bok
            self.N.append([(1 - ksi) / 2, 0, 0, (1 + ksi) / 2])  #lewy bok
               
def calcHbc(surface,element_nodes,alfa,npc):
    gaussIntegration=GaussIntegration(npc)
    Hbc=[[0.0 for _ in range(4)] for _ in range(4)]   
    
    for side in range(4): #po kazdej scianie
        for point in range(npc): #po punktach calkowania
            N=surface.N[side] 
            x1,y1=element_nodes[side].x,element_nodes[side].y
            x2, y2 = element_nodes[(side+ 1) % 4].x, element_nodes[(side + 1) % 4].y
            length = math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)
            
            #dlugosc sciany
            detJ = length / 2.0  # Wyznacznik jakobianu dla 1D
            for point, weight in enumerate(gaussIntegration.w):
                N=surface.N[side]
            

            for i in range(4):
                for j in range(4):
                    Hbc[i][j] += alfa * N[i] * N[j] * weight * detJ
      
    print("Hbc")              
    for row in Hbc:
        print(row)

    return Hbc

       


#Grid dla przykladu z prezentacji do jakobianow---------------------------
grid_przyklad=Grid(4,1) #4 nodes 1 element 
# grid_przyklad.nodes.append(Node(0, 0))          
# grid_przyklad.nodes.append(Node(0.025, 0))      
# grid_przyklad.nodes.append(Node(0.025, 0.025))  
# grid_przyklad.nodes.append(Node(0, 0.025))

# grid_przyklad.nodes.append(Node(0.01, -0.01))          
# grid_przyklad.nodes.append(Node(0.025, 0))      
# grid_przyklad.nodes.append(Node(0.025, 0.025))  
# grid_przyklad.nodes.append(Node(0, 0.025))


grid_przyklad.elements.append(Element([1, 2, 3, 4]))

grid_przyklad.display_nodes()
grid_przyklad.display_elements()
           



#---LAB 3 Jakobian---------------------------------------------

npc=9
weight=1.0 
elem_univ=ElementUniv(npc)

#Jakobian calc for every element
#DLA PRZYKLADOW TESTOWYCH----------------------------------------------------------
# for element in grid_przyklad.elements:
#     element_nodes = [grid_przyklad.nodes[id-1] for id in element.ID]  # ID węzłów zaczyna się od 1
#     jakobians = []
   
#     #Jakobian calc for every point of integration
#     for i in range(npc):
#         jakobian = Jakobian(element_nodes, elem_univ, npc)
#         jakobians.append(jakobian) #add to jakobians[]


    
#     #element stores our computed jakobian
#     element.Jakobian=jakobians
#     weights = [(1, 1), (1, 1), (1, 1), (1, 1)]
#     for j, jakobian in enumerate(jakobians):
#         print(f"Jakobian dla elementu: {element.ID}:  \nw punkcie całkowania pc{j+1}: ")
#         for row in jakobian.J:
#             print(row)
#         print("\nInverted Jakobian: ")
#         for row in jakobian.J1:
#             print(row)
#         print(f"detJ: ")
#         print(jakobian.detJ)
#         print("\n")
        
        # Wywołanie calculate_H tylko dla przykładu

            
            

         


file = open('Test1_4_4.txt', 'r')
# file = open('Test2_4_4_MixGrid.txt', 'r')
#file=open('Test3_31_31_kwadrat.txt','r')
lines=file.readlines()
file.close()

data=GlobalData(lines)
grid = Grid(data.nN,data.nE)#tworzenie obiektu grid        

load_data_from_file(lines,grid) 
k=25#conduvtivity dla pliku Test1_4_4 txt    
alfa=300 #wsp wymiany ciepla
   
# grid.display_nodes()
# grid.display_elements()

# Obliczenia jakobianu dla elementów wczytanych z pliku
# Obliczenia jakobianów dla elementów wczytanych z pliku

#global H

global_H=[[0 for _ in range(grid.nN)] for _ in range(grid.nN)]

surface = Surface(npc) #powierzchnia

for element in grid.elements:
    element_nodes = [grid.nodes[id-1] for id in element.ID]  # ID węzłów zaczyna się od 1
    jakobians = []
    
    # Obliczanie jakobianów dla każdego punktu całkowania
    for i in range(npc):
        jakobian = Jakobian(element_nodes, elem_univ, npc)
        jakobians.append(jakobian) 
    # Przechowywanie jakobianow w elemencie
    element.Jakobian = jakobians
    
  
    print(f"Obliczanie macierzy H dla elementu {element.ID}")   
    local_H=calcH(jakobian, elem_univ, k, npc) 
    Hbc=calcHbc(surface,element_nodes,alfa,npc)
    
    local_H+=Hbc
    
    agregation(global_H,element,local_H)
    
    
    
    
    
    print("Global H matrix: ")
    for row in global_H:
        print(row)
    
    
    
    
  
#     # Wyświetlanie wyników dla elementów wczytanych z pliku
#     for j, jakobian in enumerate(jakobians):
#         print(f"Jakobian dla elementu: {element.ID}: \nw punkcie całkowania pc{j+1}: ")
#         for row in jakobian.J:
#             print(row)
#         print("\nInverted Jakobian: ")
#         for row in jakobian.J1:
#             print(row)
#         print(f"detJ: {jakobian.detJ}\n")


#------DLA PRZYKLADOW TESTOWYCH----------------------------------------------
# for jakobian in element.Jakobian:  # Loop through each Jacobian
#         calcH(jakobian, elem_univ, k, npc) 


