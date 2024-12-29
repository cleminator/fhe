import poly
import util
import copy

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
    
    def cadd(self, other):
        pass
    
    def cmult(self, other):
        pass
    
    def __add__(self, other):
        b = self.b + other.b
        a = self.a + self.a
        return Ciphertext(b, a, self.q0, self.delta, self.l)
    
    def __sub__(self, other):
        b = self.b - other.b
        a = self.a - other.a
        return Ciphertext(b, a, self.q0, self.delta, self.l)
       
    def __mult__(self, other):
        # Other needs to be a tuple of (Ciphertext, evk)
        # Example: cmult = c1 * (c2, evk)
        #print("MULT")
        c1 = copy.deepcopy(self)
        c2 = copy.deepcopy(other[0])
        evk = copy.deepcopy(other[1])
        
        # Check if levels are different between both ciphertexts
        if c1.l > c2.l:
            c1.mod_reduction(self.delta)            
        elif self.l < other[0].l:
            c2.mod_reduction(self.delta)
        
        d0 = c1.b * c2.b
        d1 = c1.a * c2.b + c2.a * c1.b
        d2 = c1.a * c2.a
        
        cmult0 = d0 + round((1/self.P) * d2 * evk[0])
        cmult1 = d1 + round((1/self.P) * d2 * evk[1])
        
        cmult = Ciphertext(cmult0, cmult1, self.q0, self.delta, self.l)
        cmult.rescale()
        return cmult
        
        # Perform multiplication with relinearization
        # Rescale
        # Return
    