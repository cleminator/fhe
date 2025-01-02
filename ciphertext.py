from poly import Polynomial
import util
import copy
import ckks

class Ciphertext:
    def __init__(self, b, a, P, q0, delta, L):
        self.b = b
        self.a = a
        
        self.P = P
        
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
    
    def add_constant(self, other):
        # self: ciphertext; other: constant
        b = self.b + other
        a = self.a
        return Ciphertext(b, a, self.P, self.q0, self.delta, self.l)
        
    def add_ciph(self, other):
        b = self.b + other.b
        a = self.a + self.a
        return Ciphertext(b, a, self.P, self.q0, self.delta, self.l)
    
    def __add__(self, other):
        if isinstance(other, Ciphertext):
            return self.add_ciph(other)
        elif isinstance(other, Polynomial):
            return self.add_constant(other)
        
    
    def sub_constant(self, other):
        b = self.b - other
        a = self.a
        return Ciphertext(b, a, self.P, self.q0, self.delta, self.l)
        
    def sub_ciph(self, other):
        b = self.b - other.b
        a = self.a - other.a
        return Ciphertext(b, a, self.P, self.q0, self.delta, self.l)
    
    def __sub__(self, other):
        if isinstance(other, Ciphertext):
            return self.sub_ciph(other)
        elif isinstance(other, Polynomial):
            return self.sub_constant(other)
       
    def mult_constant(self, other):
        b = self.b * other
        a = self.a * other
        cmult = Ciphertext(b, a, self.P, self.q0, self.delta, self.l)
        cmult.rescale()
        return cmult
    
    def mult_ciph(self, other):
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
        
        cmult0 = d2 * evk[0] * (1/self.P)
        print(type(1/self.P), (1/self.P))
        cmult0.coeffs = [round(c0) for c0 in cmult0.coeffs]
        cmult0 += d0
        #cmult1 = d1 + round(d2 * evk[1] * (1/self.P))
        cmult1 = d2 * evk[1] * (1/self.P)
        cmult1.coeffs = [round(c1) for c1 in cmult1.coeffs]
        cmult1 += d1
        
        cmult = Ciphertext(cmult0, cmult1, self.P, self.q0, self.delta, self.l)
        cmult.rescale()
        return cmult
    
    def __mul__(self, other):
        # First option: Ciphertext-Ciphertext Multiplication
        #     Other needs to be a tuple of (Ciphertext, evk)
        #     Example: cmult = c1 * [c2, evk]
        # Second option: Ciphertext-Constant Mult
        #     Other needs to be an int or float
        if isinstance(other, list):
            if len(other) == 2 and isinstance(other[0], Ciphertext):
                return self.mult_ciph(other)
        elif isinstance(other, Polynomial):
            return self.mult_constant(other)
        else:
            return NotImplemented
        
    