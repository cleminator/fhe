import random
import ntt
import util

class Polynomial:
    """
    This class represents a polynomial of the integer polynomial ring Zq[X]/(X^N + 1)
    The coefficients are integers modulo q, and are represented as a list with the highest order coefficient at index 0
    """
        
    def __init__(self, coeffs, q0, delta, L):
        self.q0 = q0
        self.delta = delta
        self.L = L
        #self.q = q0 * delta**L
        self.n = len(coeffs)
        self.coeffs = [self.mod(c) for c in coeffs]
        #print(q, self.q)
        #if (2*self.n) % (q-1) != 0:
        #    raise "2n must divide q-1 so that a 2nth primitive root of unity exists!"
    
    def ql(self):
        return self.q0 * self.delta**self.L
    
    def __str__(self):
        """Return a string representation of the polynomial."""
        
        return "".join(str(self.coeffs)) + " (q = " + str(self.q0) + " * " + str(self.delta) + "^" + str(self.L) + " = " + str(self.ql()) + ")"
    
    def mod(self, val):
        z_mod_q = val % self.ql()
        if z_mod_q > self.ql() / 2:
            z_mod_q -= self.ql()
        return z_mod_q
    
    def __add__(self, other):
        max_len = max(self.n, other.n)
        result = [0] * max_len
        
        if self.ql() != other.ql():
            print("Modulus q must be the same for both polynomials")
            return None
        
        for i in range(0, max_len):
            c1 = self.coeffs[i] if i < self.n else 0
            c2 = other.coeffs[i] if i < other.n else 0
            result[i] = self.mod(c1 + c2)# % self.q
        return Polynomial(result, self.q0, self.delta, self.L)
    
    def __sub__(self, other):
        max_len = max(self.n, other.n)
        result = [0] * max_len
        
        if self.ql() != other.ql():
            print("Modulus q must be the same for both polynomials")
            return None
        
        for i in range(0, max_len):
            c1 = self.coeffs[i] if i < self.n else 0
            c2 = other.coeffs[i] if i < other.n else 0
            result[i] = self.mod(c1 - c2)# % self.q
        return Polynomial(result, self.q0, self.delta, self.L)
    
    def scalar_mult(self, scalar):
        #res = self
        cfs = [scalar * c for c in self.coeffs]
        
        return Polynomial(cfs, self.q0, self.delta, self.L)
    
    def __mul__(self, other):
        """
        a = ntt.ntt_psi(self)
        b = ntt.ntt_psi(other)
         
        coeffs = []
        for i in range(0, a.n):
            coeffs.append(self.mod(a.coeffs[i] * b.coeffs[i]))# % a.q)
        c = Polynomial(coeffs, a.q)
        res = ntt.intt_psi(c)
        return res
        """
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

        return Polynomial(result, self.q0, self.delta, self.L)
    
    def mod_exp(self, base, exp):
        result = 1
        while exp > 0:
            if exp % 2 == 1:  # If exp is odd
                result = self.mod(result * base)# % self.q
            base = self.mod(base * base)# % self.q
            exp //= 2
        return result
    
    def find_primitive_nth_root_of_unity(self, n):
        while True:
            x = random.randint(1, self.ql() - 1)
            g = self.mod_exp(x, (self.ql() - 1) / n)
            if self.mod_exp(g, n / 2) != 1:
                return g
    

    def find_2nth_root_of_unity(self, n):
        #print("2nth root")
        omega = self.find_primitive_nth_root_of_unity(n)
        #print(omega)
        for i in range(1, self.ql()):
            if self.mod_exp(i, 2) == omega:
                if self.mod_exp(i, n) == -1 or self.mod_exp(i, n) == self.ql() - 1:
                    return i
        raise "No 2nth root found"
