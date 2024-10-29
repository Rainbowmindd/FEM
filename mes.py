
class Node():
    def __init__(self,x,y):
        self.x=x
        self.y=y
        
class Element():
    def __init__(self,ID,Jakobian=None): #domyslna wartosc dla jakobianu
        
        self.ID=ID
        #lab3:
        self.Jakobian=Jakobian
        
class Grid():
    def __init__(self,nN,nE):
        self.nN=nN #number of nodes
        self.nE=nE #number of elements
        self.nodes=[]
        self.elements=[]
        
    def display_nodes(self):
        print("\nWspolrzedne wezlow: \n")
        for i, node in enumerate(self.nodes):
            print(f"Wezel {i+1}: ({node.x}, {node.y})")
            
    def display_elements(self):
        print("\nId wezlow poszczegolnych elementow: \n")
        for i, element in enumerate(self.elements):
            print(f"Element {i+1}: Wezly {element.ID}")
            
class Jakobian:
    def __init__(self):
        self.J=[[0,0],[0,0]] #2x2 matrix J
        self.J1=[[0,0],[0,0]] #2x2 matrix J^-1
        self.detJ=0 #wyznacznik
    
    def calc_jakobian(self,dN_ksi,dN_eta,element_nodes):
        x=[node.x for node in element_nodes]
        y=[node.y for node in element_nodes]
        
        #derivatives, calculate J
        #J[2][2]
        self.J[0][0] = sum(dN_ksi[i] * x[i] for i in range(npc))  # dx/dksi
        self.J[0][1] = sum(dN_ksi[i] * y[i] for i in range(npc))  # dy/dksi
        self.J[1][0] = sum(dN_eta[i] * x[i] for i in range(npc))  # dx/deta
        self.J[1][1] = sum(dN_eta[i] * y[i] for i in range(npc))  # dy/deta
        
        #wyznacznik
        self.detJ=self.J[0][0] * self.J[1][1] -self.J[0][1] *self.J[1][0]
       
        #J1[2][2]
        #inversion - odwracanie macierzy
        self.J1[0][0] = self.J[1][1] / self.detJ
        self.J1[0][1] = -self.J[0][1] / self.detJ
        self.J1[1][0] = -self.J[1][0] / self.detJ
        self.J1[1][1] = self.J[0][0] / self.detJ
       

class ElementUniv: #element uniwersalny,niezalezny od siatki
    def __init__(self,npc):
        #dN_dksi and dN_deta (pochodne funkcji ksztaltu)
        self.dN_dksi=[[0] * 4 for _ in range(npc)] #npc x 4(nodes)
        self.dN_deta=[[0] * 4 for _ in range(npc)] 
        
        #Integration points
        coordinate=1.0/(3**0.5) 
        points=[
            (-coordinate, -coordinate),
            (coordinate, -coordinate),
            (-coordinate, coordinate),
            (coordinate, coordinate)
        ]
        
        #initialize derivatives for every integration point

        for i, (ksi,eta) in enumerate(points):
            self.dN_dksi[i][0]=-0.25*(1-eta) #dN1/dksi
            self.dN_dksi[i][1]=0.25*(1-eta) #...
            self.dN_dksi[i][2]=0.25*(1+eta)
            self.dN_dksi[i][3]=-0.25*(1+eta) #dN4/dksi
            
            self.dN_deta[i][0]=-0.25*(1-ksi) #dN1/deta
            self.dN_deta[i][1]=-0.25*(1+ksi) 
            self.dN_deta[i][2]=0.25*(1+ksi) 
            self.dN_deta[i][3]=0.25*(1-ksi) #dN4/deta

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
                
            #lab 3 :
            
#LAB 4 matrix H----------------------------
def calculate_H(jakobian, elem_univ, weight, detJ):
    H = [[0.0 for _ in range(4)] for _ in range(4)] #H[4][4]
    
    #osobne partial matrices H dN_dx*dN_dx i dN_dy*dN_dy for each integration point
    HpcX=[]
    HpcY=[]
    # wyniki dla dN/dx i dN/dy
    dN_dx_table = [[0.0] * 4 for _ in range(len(elem_univ.dN_dksi))]
    dN_dy_table = [[0.0] * 4 for _ in range(len(elem_univ.dN_dksi))]

    # H w punkcie całkowania
    for i in range(len(elem_univ.dN_dksi)):
        dN_dksi = elem_univ.dN_dksi[i]
        dN_deta = elem_univ.dN_deta[i]
       

        # calc dN/dx  dN/dy
        for j in range(4):
            dN_dx = jakobian.J1[0][0] * dN_dksi[j] + jakobian.J1[0][1] * dN_deta[j]
            dN_dy = jakobian.J1[1][0] * dN_dksi[j] + jakobian.J1[1][1] * dN_deta[j]
            
            # save resu;t
            dN_dx_table[i][j] = dN_dx
            dN_dy_table[i][j] = dN_dy
        
        localX=[] #lokalne macierze dla obecnego punktu calkowania in range npc
        localY=[]

        for x in range(npc): #uwzgledniamy transpozycje
            rowX=[]
            rowY=[]
            for y in range(npc):
                rowX.append(dN_dx_table[i][x]*dN_dx_table[i][y])
                rowY.append(dN_dy_table[i][x]*dN_dy_table[i][y])
            localX.append(rowX)
            localY.append(rowY)

        HpcX.append(localX)
        HpcY.append(localY)

        # Sum H
        for j in range(4):
            for k in range(4):
                H[j][k] += (HpcX[i][j][k] + HpcY[i][j][k]) * weight * detJ *30

    # Wyświetlanie HpcX i HpcY dla każdego punktu całkowania
    for idx, (matrixX, matrixY) in enumerate(zip(HpcX, HpcY), start=1):
        print(f"\nMacierz HpcX dla punktu całkowania pc{idx}:")
        for row in matrixX:
            print(" ".join(f"{value:8.4f}" for value in row))

        print(f"\nMacierz HpcY dla punktu całkowania pc{idx}:")
        for row in matrixY:
            print(" ".join(f"{value:8.4f}" for value in row))  

    print(f"\nMacierz HpcY dla punktu całkowania pc{idx}:")
    for row in matrixY:
        print(" ".join(f"{value:8.4f}" for value in row))
    print("\nTabela dN/dx:")
    print("pc   dN1/dx   dN2/dx   dN3/dx   dN4/dx")
    for i in range(len(dN_dx_table)):
        print(f"pc{i + 1} {dN_dx_table[i][0]:8.4f} {dN_dx_table[i][1]:8.4f} {dN_dx_table[i][2]:8.4f} {dN_dx_table[i][3]:8.4f}")

    print("\nTabela dN/dy:")
    print("pc   dN1/dy   dN2/dy   dN3/dy   dN4/dy")
    for i in range(len(dN_dy_table)):
        print(f"pc{i + 1} {dN_dy_table[i][0]:8.4f} {dN_dy_table[i][1]:8.4f} {dN_dy_table[i][2]:8.4f} {dN_dy_table[i][3]:8.4f}") 
 
    print("\nMACIERZ H: ") 
    for row in H:
        print(" ".join(f"{value:8.4f}" for value in row))
        

    return H

        


#Grid dla przykladu z prezentacji do jakobianow---------------------------
grid_przyklad=Grid(4,1) #4 nodes 1 element 
grid_przyklad.nodes.append(Node(0, 0))          
grid_przyklad.nodes.append(Node(0.025, 0))      
grid_przyklad.nodes.append(Node(0.025, 0.025))  
grid_przyklad.nodes.append(Node(0, 0.025))

grid_przyklad.elements.append(Element([1, 2, 3, 4]))

grid_przyklad.display_nodes()
grid_przyklad.display_elements()
           



#---LAB 3 Jakobian---------------------------------------------
#4 integration points
npc=4
weight=1.0 
elem_univ=ElementUniv(npc)

#Jakobian calc for every element
for element in grid_przyklad.elements:
    element_nodes = [grid_przyklad.nodes[id-1] for id in element.ID]  # ID węzłów zaczyna się od 1
    jakobians = []
   
    #Jakobian calc for every point of integration
    for i in range(npc):
        jakobian=Jakobian()
        jakobian.calc_jakobian(elem_univ.dN_dksi[i], elem_univ.dN_deta[i], element_nodes)
        jakobians.append(jakobian) #add to jakobians[]


    
    #element stores our computed jakobian
    element.Jakobian=jakobian
    weights = [(1, 1), (1, 1), (1, 1), (1, 1)]
    for j, jakobian in enumerate(jakobians):
        print(f"Jakobian dla elementu: {element.ID}:  \nw punkcie całkowania pc{j+1}: ")
        for row in jakobian.J:
            print(row)
        print("\nInverted Jakobian: ")
        for row in jakobian.J1:
            print(row)
        print(f"detJ: ")
        print(jakobian.detJ)
        print("\n")
        
        # Wywołanie calculate_H tylko dla przykładu
        if element.ID == [1, 2, 3, 4]:  # Sprawdzenie, czy to ten konkretny element
            H_matrix = calculate_H(jakobian, elem_univ, weight, jakobian.detJ)
            
            

         


file = open('Test1_4_4.txt', 'r')
lines=file.readlines()
file.close()

data=GlobalData(lines)
grid = Grid(data.nN,data.nE)#tworzenie obiektu grid        

load_data_from_file(lines,grid)      
    
# grid.display_nodes()
# grid.display_elements()

# Obliczenia jakobianu dla elementów wczytanych z pliku
# for element in grid.elements:
#     element_nodes = [grid.nodes[id-1] for id in element.ID]  # ID węzłów zaczyna się od 1
#     jakobians = []
    
#     for i in range(npc):
#         jakobian = Jakobian()
#         jakobian.calc_jakobian(elem_univ.dN_dksi[i], elem_univ.dN_deta[i], element_nodes)
#         jakobians.append(jakobian)
    
#     # Wyświetlanie wyników dla elementów wczytanych z pliku
#     for j, jakobian in enumerate(jakobians):
#         print(f"Jakobian dla elementu: {element.ID}: \nw punkcie całkowania pc{j+1}: ")
#         for row in jakobian.J:
#             print(row)
#         print("\nInverted Jakobian: ")
#         for row in jakobian.J1:
#             print(row)
#         print(f"detJ: {jakobian.detJ}\n")
        


