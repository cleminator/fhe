import random
#import ntt
import util
import math


class Polynomial:
    """
    This class represents a polynomial of the integer polynomial ring Zq[X]/(X^N + 1)
    The coefficients are integers modulo q, and are represented as a list with the highest order coefficient at index 0
    """
        
    def __init__(self, coeffs, q=None):
        self.q = q #Optional modulus (encoded plaintexts do not have moduli, but ciphertexts do
        self.n = len(coeffs) #Number of coefficients
        if q:
            self.coeffs = [self.mod(c) for c in coeffs]
        else:
            self.coeffs = coeffs

    
    def __str__(self):
        """Return a string representation of the polynomial."""
        strng = "".join(str(self.coeffs))
        if self.q:
            strng += " (q = " + str(self.q) + ")"
        return strng
    
    
    ##################################################
    
    
    def mod(self, val):
        """Modulo operation with representation between -q/2 and q/2
        Source: https://eprint.iacr.org/2016/421.pdf Section 2.1 Basic"""
        z_mod_q = val % self.q
        if z_mod_q > self.q / 2:
            z_mod_q -= self.q
        return z_mod_q
    
    def rescale(self, ql):
        """Operation to get rid of the extra scaling factor after multiplying two encoded polynomials
        Source: https://eprint.iacr.org/2016/421.pdf Section 3.3"""
        self.coeffs = [math.floor(c / ql) for c in self.coeffs]
        self.q //= ql
    
    def mod_reduction(self, ql):
        """Operation to scale down the modulus without scaling down the coefficients; used to even moduli of two ciphertexts on different levels before multiplication
        Source: https://eprint.iacr.org/2016/421.pdf Section 3.3 "Homomorphic Operations of Ciphertexts at different levels"
        """
        self.q //= ql
    
    
    def mod_exp(self, base, exp):
        """ Simple modular exponentiation function"""
        if not self.q:
            raise Exception("No modulus defined!")
        result = 1
        while exp > 0:
            if exp % 2 == 1:  # If exp is odd
                result = self.mod(result * base)
            base = self.mod(base * base)
            exp //= 2
        return result
    
    ##################################################

    def solve(self, x):
        """Simple functino to solve the polynomial for x (used in decoding procedure); no optimizations performed"""
        result = 0
        for i, coeff in enumerate(self.coeffs):
            result += coeff * (x ** i)
        return result

    def __add__(self, other):
        """Adding coefficients of two polynomials"""

        from ciphertext import Ciphertext
        if isinstance(other, Ciphertext):
            return other + self
        
        max_len = max(self.n, other.n)
        result = [0] * max_len

        for i in range(0, max_len):
            c1 = self.coeffs[i] if i < self.n else 0
            c2 = other.coeffs[i] if i < other.n else 0
            result[i] = c1 + c2
            if self.q:
                result[i] = self.mod(result[i])# % self.q
        return Polynomial(result, self.q)
    
    def __sub__(self, other):
        """Subtracting coefficients of two polynomials"""
        
        max_len = max(self.n, other.n)
        result = [0] * max_len

        for i in range(0, max_len):
            c1 = self.coeffs[i] if i < self.n else 0
            c2 = other.coeffs[i] if i < other.n else 0
            result[i] = c1 - c2
            if self.q:
                result[i] = self.mod(result[i])# % self.q
        return Polynomial(result, self.q)
    
    def scalar_mult(self, scalar):
        """Multiplication of each coefficient with a constant"""
        cfs = [scalar * c for c in self.coeffs]
        return Polynomial(cfs, self.q)
    
    def __mul__(self, other):
        """Multiplication of two polynomials modulo (X^N +1); not optimized"""
        from ciphertext import Ciphertext
        if isinstance(other, Polynomial):
            return self.negacyclic_convolution(self, other)
        elif isinstance(other, int) or isinstance(other, float):
            return self.scalar_mult(other)
        elif isinstance(other, Ciphertext):
            return other * self
        else:
            return NotImplemented
            #raise Exception("Only Poly-Poly or Poly-int multiplication is posssible")


    def negacyclic_convolution(self, poly1, poly2):
        """
        Perform negacyclic convolution of two polynomials; Not yet FTT-based
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
    
    
    
    ##################################################
    
        
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
