
class Node():
    def __init__(self,x,y):
        self.x=x
        self.y=y
        
class Element():
    def __init__(self,ID):
        self.ID=ID
        
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
                
        
file = open('Test1_4_4.txt', 'r')
lines=file.readlines()
file.close()
data=GlobalData(lines)
grid = Grid(data.nN,data.nE)#tworzenie obiektu grid
           
load_data_from_file(lines,grid)
grid.display_nodes()
grid.display_elements()

       