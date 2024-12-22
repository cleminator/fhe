import numpy as np #Only use for basic stuff (exp(), complex numbers, etc.)
from poly import Polynomial
import util
import math

class CKKS:
    def __init__(self, m, scale, q):
        self.xi = np.exp(2 * np.pi * 1j / m)
        self.m = m
        self.create_sigma_R_basis()
        self.scale = scale
        self.q = q
    
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
        scaled_pi_z = self.scale * pi_z
        rounded_scaled_pi_z = self.sigma_R_discretization(scaled_pi_z)
        p = self.sigma_inverse(rounded_scaled_pi_z)
        
        coef = np.round(np.real(p.coef)).astype(int) #[round(x) for x in np.real(p.coef)]#
        p = Polynomial(coef, self.q)
        #print(p)
        return p
        
    
    def decode(self, p):
        pol = np.polynomial.Polynomial(p.coeffs)
        rescaled_p = pol / self.scale
        z = self.sigma(rescaled_p)
        pi_z = self.pi(z)
        return pi_z
        