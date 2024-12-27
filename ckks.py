import numpy as np #Only use for basic stuff (exp(), complex numbers, etc.)
from poly import Polynomial
import util
import math

class CKKS:
    def __init__(self, m, q0, delta, L, sec_dist):
        self.xi = np.exp(2 * np.pi * 1j / m)
        self.m = m
        self.create_sigma_R_basis()
        self.q0 = q0
        self.delta = delta
        self.L = L
        
        self.sec_dist = sec_dist #Distribution used for secret key term s (all error terms e, e0, e1, v, u automatically use discrete gaussian)
    
    def ql(self):
        return self.q0 * self.delta**self.L
    
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
        self.sigma_R_basis = np.array(self.vandermonde()).T
    
    def compute_basis_coordinates(self, z):
        """Computes the coordinates of a vector with respect to the orthogonal lattice basis."""
        output = np.array([np.real(np.vdot(z, b) / np.vdot(b,b)) for b in self.sigma_R_basis])
        return output

    def round_coordinates(self, coordinates):
        """Gives the integral rest."""
        coordinates = coordinates - np.floor(coordinates)
        return coordinates

    def coordinate_wise_random_rounding(self, coordinates):
        """Rounds coordinates randonmly."""
        r = self.round_coordinates(coordinates)
        f = np.array([np.random.choice([c, c-1], 1, p=[1-c, c]) for c in r]).reshape(-1)
        
        rounded_coordinates = coordinates - f
        rounded_coordinates = [int(coeff) for coeff in rounded_coordinates]
        return rounded_coordinates
        
    def sigma_R_discretization(self, z):
        """Projects a vector on the lattice using coordinate wise random rounding."""
        coordinates = self.compute_basis_coordinates(z)
        
        rounded_coordinates = self.coordinate_wise_random_rounding(coordinates)
        y = np.matmul(self.sigma_R_basis.T, rounded_coordinates)
        return y
    
    ##########
    
    def sigma_inverse(self, vec):
        van = self.vandermonde()
        coeffs = util.gaussian_elimination(van, vec)
        p = np.polynomial.Polynomial(coeffs)
        return p
    
    def sigma(self, p):
        outputs = []
        n = self.m // 2
        for i in range(n):
            root = self.xi ** (2 * i + 1)
            output = p(root)
            outputs.append(output)
        return outputs
    
    ###########
    
    def pi(self, z):
        n = self.m // 4
        return z[:n]
    
    def pi_inverse(self, z):
        z_conjugate = z[::-1]
        z_conjugate = [np.conjugate(x) for x in z_conjugate]
        return np.concatenate([z, z_conjugate])
    
    ###########
    
    def encode(self, vec):
        pi_z = self.pi_inverse(vec)
        scaled_pi_z = self.delta * pi_z
        rounded_scaled_pi_z = self.sigma_R_discretization(scaled_pi_z)
        p = self.sigma_inverse(rounded_scaled_pi_z)
        
        #coef = #np.round(np.real(p.coef)).astype(int)
        coef = [round(x) for x in np.real(p.coef)]#
        p = Polynomial(coef, self.q0, self.delta, self.L)
        #print(p)
        return p
        
    
    def decode(self, p):
        pol = np.polynomial.Polynomial(p.coeffs)
        rescaled_p = pol / self.delta
        z = self.sigma(rescaled_p)
        pi_z = self.pi(z)
        return pi_z
    
    #############
    
    def rescale(self, p):
        print("RESCALE: ")
        print(p.coeffs)
        p1 = p
        p1.coeffs = [math.floor(c / self.delta) for c in p.coeffs]
        print(p1.coeffs)
        p1.L -= 1
        print(p1)
        return p1
    
    ##############
    
    
    def keygen(self):
        match self.sec_dist:
            case 1:
                s_c = util.sample_uniform_coeffs(self.m//2, self.ql())
            case 2:
                s_c = util.sample_gaussian_coeffs(self.m//2)
            case 3:
                s_c = util.sample_uniform_ternary_coeffs(self.m//2)
            case _:
                s_c = util.sample_uniform_ternary_coeffs(self.m//2)
        s = Polynomial(s_c, self.q0, self.delta, self.L)
        a = Polynomial(util.sample_uniform_coeffs(self.m//2, self.ql()), self.q0, self.delta, self.L)
        e = Polynomial(util.sample_gaussian_coeffs(self.m//2), self.q0, self.delta, self.L)
        b = (a*s).scalar_mult(-1) + e
                
        self.sk = (1, s)
        self.pk = (b, a)
        
        # TODO: Generate EVK
    
    def encrypt(self, m):
        v = Polynomial(util.sample_gaussian_coeffs(self.m//2), self.q0, self.delta, self.L)
        e0 = Polynomial(util.sample_gaussian_coeffs(self.m//2), self.q0, self.delta, self.L)
        e1 = Polynomial(util.sample_gaussian_coeffs(self.m//2), self.q0, self.delta, self.L)
        
        #c = ((v, v) * self.pk) + (m + e0, e1)
        c0 = v * self.pk[0] + m + e0
        c1 = v * self.pk[1] + e1
        
        return (c0, c1)
    
    def decrypt(self, c):
        # c = (b, a)
        # m' = b + a * s
        return c[0] + c[1] * self.sk[1]
        