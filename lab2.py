import math

class Dane():
    def __init__(self, k):
        self.k = k  # Liczba punktów
        self.points = []  # Tablica z punktami całkowania
        self.weights = []  # Tablica z wagami
        
        self.set_points_and_weights()
        
    def set_points_and_weights(self):
          #dla N=1 , 2D 2 wartosci
        if self.k==1: #N0
            self.points=[(0,0)]
            self.weights=[4] 
        elif self.k==2: #N=1 
            self.points = [(-1 / math.sqrt(3), -1 / math.sqrt(3)),(1 / math.sqrt(3), 1 / math.sqrt(3))]
            self.weights=[1,1]
        elif self.k==3: #N=2
            self.points = [(-math.sqrt(3 / 5), -math.sqrt(3 / 5)),(0, 0),(math.sqrt(3 / 5), math.sqrt(3 / 5))]
            self.weights = [5 / 9, 8 / 9, 5 / 9]
            #elif self.k==4: i potem N4
            
    def integr_1D(self, func):
        result=0
        for i in range(len(self.points)): #waga * wart funkcji w pkt
            result+= self.weights[i] * func(self.points[i][0])
        return result
    
    def integr_2D(self,func):
        result=0
        for i in range(len(self.points)):
            for j in range(len(self.points)):
                result+= self.weights[i]*self.weights[j]*func(self.points[i][0],self.points[j][1])
        return result
                
    
def function1D(x):
    return 5 * x**2 + 3 * x + 6 
def function2D(x,y):
    return 5 * x**2 * y**2 + 3 * x * y + 6
#1D
dane_1d_2pkt=Dane(2)
result_1d_2pkt=dane_1d_2pkt.integr_1D(function1D)
print(f"Wynik calkowania 1d 2pkt :{result_1d_2pkt}")
dane_1d_3pkt=Dane(3)
result_1d_3pkt=dane_1d_3pkt.integr_1D(function1D)
print(f"Wynik calkowania 1d 3pkt :{result_1d_3pkt}")

#2D
dane_2d_2pkt=Dane(2)
result_2d_2pkt=dane_2d_2pkt.integr_2D(function2D)
print(f"Wynik calkowania 2d 2pkt :{result_2d_2pkt}")
dane_2d_3pkt=Dane(3)
result_2d_3pkt=dane_2d_3pkt.integr_2D(function2D)
print(f"Wynik calkowania 2d 3pkt :{result_2d_3pkt}")



            
                
            
        