from poly import Polynomial
from ciphertext import Ciphertext
import util
from math import e, pi, floor
import random


class CKKS:
    def __init__(self, m, P, q0, delta, L, sec_dist):
        self.xi = e ** (2 * pi * 1j / m) # E ** (2 * PI * 1j / m)
        self.m = m
        self.sigma_R_basis = self.create_sigma_R_basis()
        self.sigma_R_basis_T = self.transpose(self.sigma_R_basis)

        self.P = P

        self.q0 = q0
        self.delta = delta
        self.L = L

        self.sec_dist = sec_dist 
        #Distribution used for secret key term s (all error terms e, e0, e1, v, u automatically use discrete gaussian)
        # Source: https://eprint.iacr.org/2024/463.pdf
    
    def qL(self):
        return self.q0 * self.delta**self.L

    def transpose(self, matrix):
        rows = len(matrix)
        cols = len(matrix[0])

        transposed = [[] for _ in range(cols)]

        for row in matrix:
            for j, val in enumerate(row):
                transposed[j].append(val)

        return transposed

    def vdot(self, a, b):
        if len(a) != len(b):
            raise ValueError("Inputs must have the same length.")
        return sum((x.conjugate() if isinstance(x, complex) else x) * y for x, y in zip(a, b))

    def matmul(self, a, b):
        if len(a[0]) != len(b):
            raise ValueError("Incompatible matrix dimensions for multiplication")
        if not isinstance(b[0], list):
            b = [[bb] for bb in b]
        result = [[0] * len(b[0]) for _ in range(len(a))]
        for i in range(len(a)):  # Iterate over rows of A
            for j in range(len(b[0])):  # Iterate over columns of B
                for k in range(len(b)):  # Iterate over rows of B (columns of A)
                    result[i][j] += a[i][k] * b[k][j]

        return result

    def vandermonde(self):
        n = self.m // 2
        matrix = []
        for i in range(n):
            root = self.xi ** (2 * i + 1)
            row = []
            for j in range(n):
                row.append(root ** j)
            matrix.append(row)
        return matrix
    
    def create_sigma_R_basis(self):
        #return np.array(self.vandermonde()).T
        return self.transpose(self.vandermonde())

    def compute_basis_coordinates(self, z):
        """Computes the coordinates of a vector with respect to the orthogonal lattice basis."""
        output = [(self.vdot(z, b) / self.vdot(b,b)).real for b in self.sigma_R_basis]
        return output

    def round_coordinates(self, coordinates):
        """Gives the integral rest."""
        coordinates = [c - floor(c) for c in coordinates]
        return coordinates

    def cust_rand(self, c):
        rand = random.random()
        if rand < 1 - c:
            return c
        else:
            return 1 - c

    def coordinate_wise_random_rounding(self, coordinates):
        """Rounds coordinates randomly."""
        r = self.round_coordinates(coordinates)

        f = [self.cust_rand(c) for c in r]
        rounded_coordinates = [cc - ff for cc, ff in zip(coordinates, f)]
        rounded_coordinates = [int(coeff) for coeff in rounded_coordinates]
        return rounded_coordinates
        
    def sigma_R_discretization(self, z):
        """Projects a vector on the lattice using coordinate wise random rounding."""
        coordinates = self.compute_basis_coordinates(z)
        
        rounded_coordinates = self.coordinate_wise_random_rounding(coordinates)
        y = self.matmul(self.sigma_R_basis_T, rounded_coordinates)
        y = [yy[0] for yy in y]
        return y
    
    ##########
    
    def sigma_inverse(self, vec):
        van = self.vandermonde()
        coeffs = util.gaussian_elimination(van, vec)
        p = Polynomial(coeffs)
        return p
    
    def sigma(self, p):
        outputs = []
        n = self.m // 2
        for i in range(n):
            root = self.xi ** (2 * i + 1)
            output = p.solve(root)
            outputs.append(output)
        return outputs
    
    ###########
    
    def pi(self, z):
        n = self.m // 4
        return z[:n]
    
    def pi_inverse(self, z):
        z_conjugate = z[::-1]
        z_conjugate = [x.conjugate() for x in z_conjugate]
        return z + z_conjugate #np.concatenate([z, z_conjugate])
    
    
    ############
    
    
    def encode(self, vec):
        #Convert list of real numbers to complex numbers
        vec = [complex(vv, 0) for vv in vec]
        pi_z = self.pi_inverse(vec)
        scaled_pi_z = [self.delta * z for z in pi_z] #self.delta * pi_z
        rounded_scaled_pi_z = self.sigma_R_discretization(scaled_pi_z)
        p = self.sigma_inverse(rounded_scaled_pi_z)
        
        #coef = #np.round(np.real(p.coef)).astype(int)
        coef = [round(x.real) for x in p.coeffs]
        p = Polynomial(coef, self.qL())
        return p
        
    
    def decode(self, p):
        rescaled_p = Polynomial([c / self.delta for c in p.coeffs])
        z = self.sigma(rescaled_p)
        pi_z = self.pi(z)
        #Extract real parts of complex vector
        pi_z = [c.real for c in pi_z]
        return pi_z
    
    
    ######################################################
       
    # Function to generate public key, private key
    # Source: https://eprint.iacr.org/2016/421.pdf (Section 3.4)
    def keygen(self):
        match self.sec_dist:
            case 1:
                s_c = util.sample_uniform_coeffs(self.m//2, self.qL())
            case 2:
                s_c = util.sample_gaussian_coeffs(self.m//2)
            case 3:
                s_c = util.sample_uniform_ternary_coeffs(self.m//2)
            case _:
                s_c = util.sample_uniform_ternary_coeffs(self.m//2)
        s = Polynomial(s_c, self.qL())
        a = Polynomial(util.sample_uniform_coeffs(self.m//2, self.qL()), self.qL())
        e = Polynomial(util.sample_gaussian_coeffs(self.m//2), self.qL())
        b = a*s
        b = b * -1
        b = b + e
        
        sk = (1, s)
        pk = (b, a)
        
        return (pk, sk)
    
    # Function to generate evaluation key for Ciphertext-Ciphertext multiplication
    # Source: https://eprint.iacr.org/2016/421.pdf (Section 3.4)
    def evkeygen(self, sk):
        ap = Polynomial(util.sample_uniform_coeffs(self.m//2, self.P * self.qL()), self.P * self.qL())
        ep = Polynomial(util.sample_gaussian_coeffs(self.m//2), self.P * self.qL())

        bp = ap * sk[1]
        bp = bp * -1
        bp = bp + ep
        p_term = (sk[1] * sk[1] * self.P)
        bp = bp + p_term
        evk = (bp, ap)
        return evk
    
    def encrypt(self, m, pk):
        #v = Polynomial(util.sample_gaussian_coeffs(self.m//2), self.qL())
        v_c = util.sample_gaussian_coeffs(self.m//2)#[abs(x) for x in util.sample_uniform_ternary_coeffs(self.m//2)]
        #print(v_c)
        v = Polynomial([1]*(self.m//2), self.qL())
        e0 = Polynomial(util.sample_gaussian_coeffs(self.m//2), self.qL())
        e1 = Polynomial(util.sample_gaussian_coeffs(self.m//2), self.qL())
        
        #c = ((v, v) * self.pk) + (m + e0, e1)
        c0 = (v * pk[0]) + m + e0
        c1 = (v * pk[1]) + e1
        
        c = Ciphertext(c0, c1, self.P, self.q0, self.delta, self.L)
        return c
    
    def decrypt(self, c, sk):
        # c = (b, a)
        # m' = b + a * s
        return c.b + (c.a * sk[1])
