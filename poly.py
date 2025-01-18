
import ntt
import util
import math
from copy import copy


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
        return util.mod(val, self.q)
    
    def rescale(self, ql):
        """Operation to get rid of the extra scaling factor after multiplying two encoded polynomials
        Source: https://eprint.iacr.org/2016/421.pdf Section 3.3"""
        self.coeffs = [(c // ql) for c in self.coeffs]
        self.q //= ql
    
    def mod_reduction(self, ql):
        """Operation to scale down the modulus without scaling down the coefficients; used to even moduli of two ciphertexts on different levels before multiplication
        Source: https://eprint.iacr.org/2016/421.pdf Section 3.3 "Homomorphic Operations of Ciphertexts at different levels"
        """
        self.q //= ql
    
    

    
    ##################################################

    def solve(self, x):
        """Simple function to solve the polynomial for x (used in decoding procedure); no optimizations performed"""
        result = 0
        for i, coeff in enumerate(self.coeffs):
            result += coeff * (x ** i)
        return result

    def add(self, other):
        max_len = max(self.n, other.n)
        result = [0] * max_len

        for i in range(0, max_len):
            c1 = self.coeffs[i] if i < self.n else 0
            c2 = other.coeffs[i] if i < other.n else 0
            result[i] = c1 + c2
            if self.q:
                result[i] = self.mod(result[i])  # % self.q
        return Polynomial(result, self.q)

    def __add__(self, other):
        """Adding coefficients of two polynomials"""

        from ciphertext import Ciphertext
        if isinstance(other, Ciphertext):
            return other + self
        
        return self.add(other)

    def sub(self, other):
        max_len = max(self.n, other.n)
        result = [0] * max_len

        for i in range(0, max_len):
            c1 = self.coeffs[i] if i < self.n else 0
            c2 = other.coeffs[i] if i < other.n else 0
            result[i] = c1 - c2
            if self.q:
                result[i] = self.mod(result[i])  # % self.q
        return Polynomial(result, self.q)

    def __sub__(self, other):
        """Subtracting coefficients of two polynomials"""
        return self.sub(other)
    
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
            #raise Exception("Only Poly-Poly or Poly-int multiplication is possible")


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
    
        



class RNSLimb:
    def __init__(self, coeffs, q, root):
        self.coeffs = coeffs
        self.q = q
        self.root = root
    def __copy__(self):
        return RNSLimb(self.coeffs[:], self.q, self.root)
    def __getitem__(self, item):
        return self.coeffs[item]
    def __setitem__(self, key, value):
        self.coeffs[key] = value
    def __iter__(self):
        for c in self.coeffs:
            yield c
    def __len__(self):
        return len(self.coeffs)
    def __str__(self):
        strng = "".join(str(self.coeffs))
        strng += " (Modulus: " + str(self.q) + ", Root: " + str(self.root) + ")"
        return strng



class RNSPolynomial(Polynomial):
    # B: RNS base of P [p0, p1, ..., pk-1] (len(B) = k)
    # C: RNS base of Q [q0, q1, ..., qL] (len(C) = L+1)
    # limbs: list of RNS limbs, each limb i is a list of coeffs mod q_i

    def __init__(self, B, C, roots, is_ntt=False, coeffs=None):
        self.B = B[:]
        self.C = C[:]
        self.roots = roots
        self.limbs = []

        self.n = 0
        self.ntt_domain = is_ntt
        if coeffs:
            self.create_limbs(coeffs[:])


    def __str__(self):
        strng = "".join(str(self.get_coeffs()))
        strng += " (Modulus = " + str(self.get_modulus()) + ")"
        if self.ntt_domain:
            strng += " (NTT)"
        return strng


    def get_coeffs(self):

        if len(self.C) == 1:
            return self.limbs[0]
        if len(self.C) == 0:
            raise Exception("No limbs exist")

        # x = l1*m2*n2 + l2*m1*n1

        coeffs = [0]*len(self.limbs[0])
        l1 = self.limbs[0]
        q1 = self.C[0]

        for i in range(1, len(self.C)):
            _, m1, m2 = util.extGCD(q1, self.C[i])
            for j in range(len(l1)):
                coeffs[j] = util.mod(l1[j] * m2 * self.C[i] + self.limbs[i][j] * m1 * q1, self.get_modulus())
            l1 = coeffs
            q1 = q1 * self.C[i]
        return coeffs

    def get_modulus(self):
        q = 1
        for l in self.limbs:#range(len(self.C)):
            q *= l.q#self.C[i]
        return q

    ###########################################################################

    def create_limbs(self, coeffs):
        # Generate list of limbs from base (Q) and coefficients
        self.n = len(coeffs[:])
        for i in range(len(self.C)):
            #print("Generating limb", i)
            l = RNSLimb([c % self.C[i] for c in coeffs[:]], self.C[i], self.roots[self.C[i]])
            self.limbs.append(l) #[c % self.C[i] for c in coeffs[:]])

        if not self.ntt_domain:
            self.convert_RNS_to_NTT()
            self.ntt_domain = True

    def set_limbs(self, limbs):
        self.n = len(limbs[0])
        self.limbs = [copy(l) for l in limbs[:]] #limbs[:]]
        if not self.ntt_domain:
            self.convert_RNS_to_NTT()
        return self

    ###########################################################################

    ## These functions will transform each of the limbs to and from NTT domain

    def convert_RNS_to_NTT(self):
        for l in self.limbs:
            #print("Converting to NTT limb")
            l.coeffs = ntt.ntt_psi(l.coeffs, self.n, l.q, l.root)
        self.ntt_domain = True

    def convert_NTT_to_RNS(self):
        for l in self.limbs:
            l.coeffs = ntt.intt_psi(l.coeffs, self.n, l.q, l.root)
        self.ntt_domain = False

    ###########################################################################

    def conv(self, base1, base2):
        # Convert limbs from basis1 (e.g. B) to basis2 (e.g. C)
        pass

    def mod_up(self):
        pass

    def mod_down(self):
        pass

    def solve(self, x):
        """Simple function to solve the polynomial for x (used in decoding procedure); no optimizations performed"""
        result = 0
        for i, coeff in enumerate(self.get_coeffs()):
            result += coeff * (x ** i)
        return result

    ######################################

    def rescale(self):
        for l in self.limbs:
            limb = []
            for i in range(self.n):
                lp = l[i] - self.limbs[-1][i]
                lp = lp * util.findMultInv(self.limbs[-1].q, l.q)
                lp = util.mod(lp, l.q)
                limb.append(lp)
            l.coeffs = limb
        """
        limbs_p = []
        for lj, qj in zip(self.limbs, self.C):
            #print("Limb: ", lj, qj)
            limb = []
            for i in range(self.n):
                #print("i: ", i, lj[i], self.limbs[-1][i])
                l = lj[i]-self.limbs[-1][i]
                #print(l)
                l = l * util.findMultInv(self.C[-1], qj) # // self.C[-1]
                #print(l)
                l = util.mod(l, qj)
                #print(l)
                limb.append(l)
            #print(limb)
            limbs_p.append(RNSLimb(limb, ))
        #print(self.limbs)
        #print(limbs_p)
        """
        #self.set_limbs(limbs_p)
        self.limbs.pop()
        self.C.pop()


    def mod_reduction(self):
        self.limbs.pop()
        self.C.pop()

    ######################################

    def __add__(self, other):
        from ciphertext import RNSCiphertext
        if isinstance(other, RNSCiphertext):
            return other + self

        if self.n != other.n:
            raise Exception("Polynomials need to have equal ring dimension")


        limbs = []
        # zip() will stop at the end of the shortest list.
        # If one poly is lower level, the result will therefore automatically be reduced to that level
        for la, lb, q in zip(self.limbs, other.limbs, self.C):
            lc = [0] * self.n
            for i in range(self.n):
                lc[i] = util.mod(la[i] + lb[i], q)
            limbs.append(RNSLimb(lc, la.q, la.root))#lc)
        return RNSPolynomial(self.B, self.C[:len(limbs)], self.roots, self.ntt_domain).set_limbs(limbs)

    def __sub__(self, other):
        if self.n != other.n:
            raise Exception("Polynomials need to have equal ring dimension")

        limbs = []
        # zip() will stop at the end of the shortest list.
        # If one poly is lower level, the result will therefore automatically be reduced to that level
        for la, lb, q in zip(self.limbs, other.limbs, self.C):
            lc = [0] * self.n
            for i in range(self.n):
                lc[i] = util.mod(la[i] - lb[i], q)
            limbs.append(RNSLimb(lc, la.q, la.root))#lc)
        return RNSPolynomial(self.B, self.C[:len(limbs)], self.roots, self.ntt_domain).set_limbs(limbs)

    def scalar_mult(self, scalar):
        limbs = []
        for l, q in zip(self.limbs, self.C):
            limbs.append(RNSLimb([util.mod(scalar * c, q) for c in l], l.q, l.root))
        return RNSPolynomial(self.B, self.C, self.roots, self.ntt_domain).set_limbs(limbs)

    def negacyclic_convolution(self, poly1, poly2):
        limbs = []
        for la, lb in zip(poly1.limbs, poly2.limbs):
            n = len(la)
            m = len(lb)

            if n != m:
                raise ValueError("Negacyclic convolution requires polynomials of the same length.")

            #if not poly1.ntt_domain:
            #    poly1.convert_RNS_to_NTT()
            #if not poly2.ntt_domain:
            #    poly2.convert_RNS_to_NTT()

            result = [0] * n

            for i in range(self.n):
                result[i] = la[i] * lb[i]

            limbs.append(RNSLimb([util.mod(r, la.q) for r in result], la.q, la.root))
        return RNSPolynomial(poly1.B, poly1.C, poly1.roots, poly1.ntt_domain).set_limbs(limbs)

    def __mul__(self, other):
        """Multiplication of two polynomials modulo (X^N +1); not optimized"""
        from ciphertext import RNSCiphertext
        if isinstance(other, Polynomial):
            return self.negacyclic_convolution(self, other)
        elif isinstance(other, int) or isinstance(other, float):
            return self.scalar_mult(other)
        elif isinstance(other, RNSCiphertext):
            return other * self
        else:
            return NotImplemented
