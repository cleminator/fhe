import poly
import util

class Ciphertext:
    def __init__(self, b, a, q0, delta, L):
        self.b = b
        self.a = a
        
        self.q0 = q0
        self.delta = delta
        self.l = L
    
    def __str__(self):
        return "Ciphertext (q0: " + str(self.q0) + ", delta: " + str(self.delta) + ", l: " + str(self.l) + ")"
            
    ############
    
    def rescale(self):
        self.b.rescale(self.delta)
        self.a.rescale(self.delta)
        self.l -= 1
    
    ############
    
    def __add__(self, other):
        pass
    
    def __sub__(self, other):
        pass
       
    def __mult__(self, other):
        # Implement mult
        self.rescale()
        pass
    
    