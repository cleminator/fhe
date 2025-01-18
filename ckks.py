from poly import Polynomial, RNSPolynomial
from ciphertext import Ciphertext, RNSCiphertext
import util
from math import e, pi, prod
import numpy as np
import ntt


class CKKS:
    def __init__(self, N, P, q0, delta, L, sec_dist):
        self.m = 2*N
        self.xi = e ** (2 * pi * 1j / self.m)  # Used for encoding
        self.sigma_R_basis = self.create_sigma_R_basis() # used for encoding
        self.sigma_R_basis_T = util.transpose(self.sigma_R_basis) # precompute transposition

        self.P = P #Used to scale evaluation key to avoid large error terms

        self.q0 = q0 # First mod
        self.delta = delta # Scale mod
        self.L = L # Number of available levels


        self.sec_dist = sec_dist
        """Distribution used for secret key term s (all error terms e, e0, e1, v, u automatically use discrete gaussian)
        Source: https://eprint.iacr.org/2024/463.pdf """
    
    def qL(self):
        """Return the full modulus"""
        return self.q0 * self.delta**self.L


    def vandermonde(self):
        """Computes the Vandermonde matrix from an m-th root of unity.
        Source: https://blog.openmined.org/ckks-explained-part-1-simple-encoding-and-decoding/"""
        n = self.m // 2
        matrix = []
        for i in range(n):
            root = self.xi ** (2 * i + 1)
            row = []
            for j in range(n):
                row.append(root ** j)
            matrix.append(row)
        return matrix

    ########

    def create_sigma_R_basis(self):
        """Creates the basis (sigma(1), sigma(X), ..., sigma(X** N-1))
        Source: https://blog.openmined.org/ckks-explained-part-2-ckks-encoding-and-decoding/"""
        return util.transpose(self.vandermonde())

    def compute_basis_coordinates(self, z):
        """Computes the coordinates of a vector with respect to the orthogonal lattice basis
        Source: https://blog.openmined.org/ckks-explained-part-2-ckks-encoding-and-decoding/"""
        output = [(util.vdot(z, b) / util.vdot(b,b)).real for b in self.sigma_R_basis]
        return output
        
    def sigma_R_discretization(self, z):
        """Projects a vector on the lattice using coordinate wise random rounding.
        Source: https://blog.openmined.org/ckks-explained-part-2-ckks-encoding-and-decoding/"""
        coordinates = self.compute_basis_coordinates(z)
        
        rounded_coordinates = util.coordinate_wise_random_rounding(coordinates)
        y = util.matmul(self.sigma_R_basis_T, rounded_coordinates)
        y = [yy[0] for yy in y] # "Flattens" the vector from the from [[1], [2], [3]] to [1, 2, 3]
        return y
    
    ##########
    
    def sigma_inverse(self, vec):
        """Encodes the vector b in a polynomial using an M-th root of unity.
        Source: https://blog.openmined.org/ckks-explained-part-1-simple-encoding-and-decoding/"""
        van = self.vandermonde()
        #coeffs = util.gaussian_elimination(van, vec)
        coeffs = np.linalg.solve(van, vec)
        p = Polynomial(coeffs)
        return p
    
    def sigma(self, p):
        """Decodes a polynomial by applying it to the M-th roots of unity.
        Source: https://blog.openmined.org/ckks-explained-part-1-simple-encoding-and-decoding/"""
        outputs = []
        n = self.m // 2
        for i in range(n):
            root = self.xi ** (2 * i + 1)
            output = p.solve(root)
            outputs.append(output)
        return outputs
    
    ###########
    
    def pi(self, z):
        """Projects a vector of H into C^{N/2}
        Source: https://blog.openmined.org/ckks-explained-part-2-ckks-encoding-and-decoding/"""
        n = self.m // 4
        return z[:n]
    
    def pi_inverse(self, z):
        """Expands a vector of C^{N/2} by expanding it with its complex conjugate
        Source: https://blog.openmined.org/ckks-explained-part-2-ckks-encoding-and-decoding/"""
        z_conjugate = z[::-1]
        z_conjugate = [x.conjugate() for x in z_conjugate]
        return z + z_conjugate
    
    
    ############
    
    
    def encode(self, vec):
        """Encodes a vector by expanding it first to H, scale it, project it on the lattice of sigma(R), and performs sigma inverse.
        Source: https://blog.openmined.org/ckks-explained-part-2-ckks-encoding-and-decoding/"""

        vec = [complex(vv, 0) for vv in vec] #Convert list of real numbers to complex numbers
        pi_z = self.pi_inverse(vec)
        scaled_pi_z = [self.delta * z for z in pi_z] #self.delta * pi_z
        rounded_scaled_pi_z = self.sigma_R_discretization(scaled_pi_z)
        p = self.sigma_inverse(rounded_scaled_pi_z)

        coef = [round(x.real) for x in p.coeffs]
        p = Polynomial(coef, self.qL())
        return p
        
    
    def decode(self, p):
        """Decodes a polynomial by removing the scale, evaluating on the roots, and project it on C^(N/2)
        Source: https://blog.openmined.org/ckks-explained-part-2-ckks-encoding-and-decoding/"""
        rescaled_p = Polynomial([c / self.delta for c in p.coeffs])
        z = self.sigma(rescaled_p)
        pi_z = self.pi(z)
        #Extract real parts of complex vector
        pi_z = [c.real for c in pi_z]
        return pi_z
    
    
    ######################################################


    def keygen(self):
        """Function to generate public key, private key
        Source: https://eprint.iacr.org/2016/421.pdf (Section 3.4)"""
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
        
        return pk, sk
    

    def evkeygen(self, sk):
        """ Function to generate evaluation key for Ciphertext-Ciphertext multiplication
        Source: https://eprint.iacr.org/2016/421.pdf (Section 3.4)"""
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
        """Encrypts a previously encoded plaintext into a ciphertext using the public key
        Source: https://eprint.iacr.org/2016/421.pdf (Section 3.4)"""
        #v = Polynomial(util.sample_gaussian_coeffs(self.m//2), self.qL())
        v = Polynomial([1]*(self.m//2), self.qL()) # This polynomial with 1s as coefficients is a placeholder until I understand the purpose of the term "v"
        e0 = Polynomial(util.sample_gaussian_coeffs(self.m//2), self.qL())
        e1 = Polynomial(util.sample_gaussian_coeffs(self.m//2), self.qL())
        
        #c = ((v, v) * self.pk) + (m + e0, e1)
        c0 = (v * pk[0]) + m + e0
        c1 = (v * pk[1]) + e1
        
        c = Ciphertext(c0, c1, self.P, self.q0, self.delta, self.L)
        return c
    
    def decrypt(self, c, sk):
        """Decrypts a ciphertext using the secret key and returns a plaintext polynomial
        Source: https://eprint.iacr.org/2016/421.pdf (Section 3.4)"""
        return c.b + (c.a * sk[1])


class RNSCKKS(CKKS):
    def __init__(self, N, p, k, q0, q, L, sec_dist):#C, q, sec_dist):
        self.B = []
        self.C = []
        self.m = 2*N
        self.xi = e ** (2 * pi * 1j / self.m)  # Used for encoding
        #print("Create sigma_R_basis")
        self.sigma_R_basis = self.create_sigma_R_basis()  # used for encoding
        self.sigma_R_basis_T = util.transpose(self.sigma_R_basis)  # precompute transposition

        self.k = k
        self.p = 2**p
        #print("Create B")
        self.B = self.generate_basis(p, k)

        self.L = L
        #print("Create C")
        self.C = self.generate_basis(q0, 1)
        self.C += self.generate_basis(q, L) #self.generate_base_C(q0, q, L)
        self.q = 2**q # All moduli q_i should be as close to this q as possible, to reduce the approximation error during rescaling
        self.q0 = 2**q0
        # q will be used for the scaling during encoding and decoding

        self.roots = {} #Will contain a map of modulus->root for each q and p
        #print("Create roots")
        self.generate_roots()
        self.sec_dist = sec_dist

    def qL(self):
        """Return the full modulus"""
        return prod(self.C)

    def generate_roots(self):
        for q in self.C:
            self.roots[q] = util.find_2nth_root_of_unity(self.m//2, q)
        for p in self.B:
            self.roots[p] = util.find_2nth_root_of_unity(self.m//2, p)

    def generate_basis(self, bitsize, length):
        b = []
        next_prime = 2**bitsize
        for i in range(length):
            while True:
                next_prime = util.find_next_prime(next_prime)
                if util.mod(next_prime, self.m) == 1 and next_prime not in self.B and next_prime not in self.C:
                    b.append(next_prime)
                    next_prime += 1
                    break
                else:
                    next_prime += 1
        return b



    def sigma_inverse(self, vec):
        """Encodes the vector b in a polynomial using an M-th root of unity.
        Source: https://blog.openmined.org/ckks-explained-part-1-simple-encoding-and-decoding/"""
        van = self.vandermonde()
        coeffs = util.gaussian_elimination(van, vec)
        #p = Polynomial(coeffs)
        #p = RNSPolynomial(self.B, self.C, [x.real for x in coeffs])
        return [c.real for c in coeffs]


    def encode(self, vec):
        """Encodes a vector by expanding it first to H, scale it, project it on the lattice of sigma(R), and performs sigma inverse.
        Source: https://blog.openmined.org/ckks-explained-part-2-ckks-encoding-and-decoding/"""
        vec = [complex(vv, 0) for vv in vec]  # Convert list of real numbers to complex numbers
        pi_z = self.pi_inverse(vec)
        scaled_pi_z = [self.q * z for z in pi_z]  # self.delta * pi_z
        rounded_scaled_pi_z = self.sigma_R_discretization(scaled_pi_z)
        p = self.sigma_inverse(rounded_scaled_pi_z)

        coef = [round(x) for x in p]
        #p = Polynomial(coef, self.qL())

        p = RNSPolynomial(self.B, self.C, self.roots, False, coef)

        return p

    def decode(self, p):
        """Decodes a polynomial by removing the scale, evaluating on the roots, and project it on C^(N/2)
        Source: https://blog.openmined.org/ckks-explained-part-2-ckks-encoding-and-decoding/"""
        #rescaled_p = Polynomial([c / self.delta for c in p.coeffs])
        pcopy = RNSPolynomial(p.B[:], p.C[:], p.roots, p.ntt_domain)
        pcopy.set_limbs(p.limbs[:])
        if pcopy.ntt_domain:
            pcopy.convert_NTT_to_RNS()
        rescaled_p = Polynomial([c / self.q for c in pcopy.get_coeffs()])
        z = self.sigma(rescaled_p)
        pi_z = self.pi(z)
        # Extract real parts of complex vector
        pi_z = [c.real for c in pi_z]
        return pi_z

    ###########################

    def keygen(self):
        """Function to generate public key, private key
        Source: https://eprint.iacr.org/2016/421.pdf (Section 3.4)"""
        match self.sec_dist:
            case 1:
                s_c = util.sample_uniform_coeffs(self.m // 2, self.qL())
            case 2:
                s_c = util.sample_gaussian_coeffs(self.m // 2)
            case 3:
                s_c = util.sample_uniform_ternary_coeffs(self.m // 2)
            case _:
                s_c = util.sample_uniform_ternary_coeffs(self.m // 2)
        #s = RNSPolynomial(self.B, self.C, self.roots, coeffs=s_c)s
        #a = RNSPolynomial(self.B, self.C, self.roots, coeffs=util.sample_uniform_coeffs(self.m // 2, self.qL()))
        #print("Sampling")
        s = RNSPolynomial(self.B, self.C, self.roots, coeffs=[1]*(self.m // 2))
        a = RNSPolynomial(self.B, self.C, self.roots, coeffs=[1]*(self.m // 2))
        e = RNSPolynomial(self.B, self.C, self.roots, coeffs=util.sample_gaussian_coeffs(self.m // 2))
        #print("Mult poly")
        b = a * s
        #print("mult const")
        b = b * -1
        #b = b + e

        sk = (1, s)
        pk = (b, a)

        return pk, sk

    def encrypt(self, m, pk):
        """Encrypts a previously encoded plaintext into a ciphertext using the public key
        Source: https://eprint.iacr.org/2016/421.pdf (Section 3.4)"""
        # v = Polynomial(util.sample_gaussian_coeffs(self.m//2), self.qL())
        v = RNSPolynomial(self.B, self.C, self.roots, coeffs=[1] * (self.m // 2))  # This polynomial with 1s as coefficients is a placeholder until I understand the purpose of the term "v"
        e0 = RNSPolynomial(self.B, self.C, self.roots, coeffs=util.sample_gaussian_coeffs(self.m // 2))
        e1 = RNSPolynomial(self.B, self.C, self.roots, coeffs=util.sample_gaussian_coeffs(self.m // 2))

        # c = ((v, v) * self.pk) + (m + e0, e1)
        c0 = (v * pk[0]) + m# + e0
        c1 = (v * pk[1])# + e1

        c = RNSCiphertext(c0, c1, self.B, self.C, self.p, self.q0, self.q, self.roots) #c0, c1, self.P, self.q0, self.delta, self.L)
        return c

    def decrypt(self, c, sk):
        """Decrypts a ciphertext using the secret key and returns a plaintext polynomial
        Source: https://eprint.iacr.org/2016/421.pdf (Section 3.4)"""
        return c.b + (c.a * sk[1])
