
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
        #J[4][4]
        self.J[0][0] = sum(dN_ksi[i][0] * x[i] for i in range(4))  # dx/dksi
        self.J[0][1] = sum(dN_ksi[i][0] * y[i] for i in range(4))  # dy/dksi
        self.J[1][0] = sum(dN_eta[i][0] * x[i] for i in range(4))  # dx/deta
        self.J[1][1] = sum(dN_eta[i][0] * y[i] for i in range(4))  # dy/deta
        
        #wyznacznik
        self.detJ=self.J[0][0] * self.J[1][1] -self.J[0][1] *self.J[1][0]
       
        #J1[4][4]
        #inversion - odwracanie macierzy
        self.J1[0][0] = self.J[1][1] / self.detJ
        self.J1[0][1] = -self.J[0][1] / self.detJ
        self.J1[1][0] = -self.J[1][0] / self.detJ
        self.J1[1][1] = self.J[0][0] / self.detJ
       

class ElementUniv: #element uniwersalny,niezalezny od siatki
    def __init__(self,npc):
        #dN_dksi and dN_deta (pochodne funkcji ksztaltu)
        self.dN_dksi=[[0] * 4 for _ in range(npc)] #npc x 4
        self.dN_deta=[[0] * 4 for _ in range(npc)] 
        
        #Integration points
        coordinate=1.0/(3**0.5) 
        points=[
            (-coordinate, -coordinate),
            (coordinate, -coordinate),
            (coordinate, coordinate),
            (-coordinate, coordinate)
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
            npc=(0,0)
                
        
file = open('Test1_4_4.txt', 'r')
lines=file.readlines()
file.close()
data=GlobalData(lines)
grid = Grid(data.nN,data.nE)#tworzenie obiektu grid
           
load_data_from_file(lines,grid)


#---LAB 3 Jakobian---------------------------------------------
#4 integration points
npc=4
elem_univ=ElementUniv(npc)

#Jakobian calc for every element
for element in grid.elements:
    
    element_nodes=[grid.nodes[id-1] for id in element.ID] #id-1 bo nodeID
    #zaczyna się od 1
    
    #Jakobian calc for every element
    jakobian=Jakobian()
    jakobian.calc_jakobian(elem_univ.dN_dksi, elem_univ.dN_deta, element_nodes)

    
    #element stores our computed jakobian
    element.Jakobian=jakobian
    
    print(f"Jakobian dla elementu: {element.ID}: ")
    for row in element.Jakobian.J:
        print(row)
    
grid.display_nodes()
grid.display_elements()