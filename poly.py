import random
import ntt
import util

class Polynomial:
    """
    This class represents a polynomial of the integer polynomial ring Zq[X]/(X^N + 1)
    The coefficients are integers modulo q, and are represented as a list with the highest order coefficient at index 0
    """
        
    def __init__(self, coeffs, q=None):
        self.q = q
        #self.delta = delta
        #self.L = L
        #self.q = q0 * delta**L
        self.n = len(coeffs)
        if q:
            self.coeffs = [self.mod(c) for c in coeffs]
        else:
            self.coeffs = coeffs
        #print(q, self.q)
        #if (2*self.n) % (q-1) != 0:
        #    raise "2n must divide q-1 so that a 2nth primitive root of unity exists!"
    
    
    
    def __str__(self):
        """Return a string representation of the polynomial."""
        strng = "".join(str(self.coeffs))
        if self.q:
            strng += " (q = " + str(self.q) + ")"
        return strng
    
    def mod(self, val):
        z_mod_q = val % self.q
        if z_mod_q > self.q / 2:
            z_mod_q -= self.q
        return z_mod_q
    
    def rescale(self, ql):
        #print("RESCALE: ")
        #print(p.coeffs)
        #p1 = p
        self.coeffs = [math.floor(c / q) for c in self.coeffs]
        #print(p1.coeffs)
        #p1.L -= 1
        self.q //= ql
        #print(p1)
        #return p1
    
    def __add__(self, other):
        max_len = max(self.n, other.n)
        result = [0] * max_len
        
        if self.q != other.q:
            print("Modulus q must be the same for both polynomials")
            return None
        
        for i in range(0, max_len):
            c1 = self.coeffs[i] if i < self.n else 0
            c2 = other.coeffs[i] if i < other.n else 0
            result[i] = c1 + c2
            if self.q:
                result[i] = self.mod(result[i])# % self.q
        return Polynomial(result, self.q)
    
    def __sub__(self, other):
        max_len = max(self.n, other.n)
        result = [0] * max_len
        
        if self.q != other.q:
            print("Modulus q must be the same for both polynomials")
            return None
        
        for i in range(0, max_len):
            c1 = self.coeffs[i] if i < self.n else 0
            c2 = other.coeffs[i] if i < other.n else 0
            result[i] = c1 - c2
            if self.q:
                result[i] = self.mod(result[i])# % self.q
        return Polynomial(result, self.q)
    
    def scalar_mult(self, scalar):
        #res = self
        cfs = [scalar * c for c in self.coeffs]
        
        return Polynomial(cfs, self.q)
    
    def __mul__(self, other):
        return self.negacyclic_convolution(self, other)

    def negacyclic_convolution(self, poly1, poly2):
        """
        Perform negacyclic convolution of two polynomials.
        :param poly1: Polynomial object representing the first polynomial.
        :param poly2: Polynomial object representing the second polynomial.
        :return: Polynomial object representing the result of the negacyclic convolution.
        """
        n = len(poly1.coeffs)
        m = len(poly2.coeffs)

        if n != m:
            raise ValueError("Negacyclic convolution requires polynomials of the same length.")

        # Initialize the result coefficients
        result = [0] * n

        # Perform direct computation of negacyclic convolution
        for i in range(n):
            for j in range(n):
                index = (i + j) % n
                value = poly1.coeffs[i] * poly2.coeffs[j]
                if i + j >= n:
                    result[index] -= value  # Negacyclic wraparound
                else:
                    result[index] += value

        return Polynomial(result, self.q)
    
    def mod_exp(self, base, exp):
        if not self.q:
            raise Exception("No modulus defined!")
        result = 1
        while exp > 0:
            if exp % 2 == 1:  # If exp is odd
                result = self.mod(result * base)# % self.q
            base = self.mod(base * base)# % self.q
            exp //= 2
        return result
    
    def find_primitive_nth_root_of_unity(self, n):
        while True:
            x = random.randint(1, self.q - 1)
            g = self.mod_exp(x, (self.q - 1) / n)
            if self.mod_exp(g, n / 2) != 1:
                return g
    

    def find_2nth_root_of_unity(self, n):
        #print("2nth root")
        omega = self.find_primitive_nth_root_of_unity(n)
        #print(omega)
        for i in range(1, self.q):
            if self.mod_exp(i, 2) == omega:
                if self.mod_exp(i, n) == -1 or self.mod_exp(i, n) == self.q - 1:
                    return i
        raise "No 2nth root found"
