from poly import Polynomial, RNSPolynomial
import ckks
from ntt import ntt_psi, intt_psi
import util
from timeit import default_timer as timer

N = 2**2
q = 30
q0 = 40
L = 4
k = 1
#P = 2**90
p = 18

"""
rnsckks = ckks.RNSCKKS(N, p, k, q0, q, L, 2)

B = []
C = rnsckks.generate_basis(20, 2)
print(C)
roots = {}
roots[C[0]] = util.find_2nth_root_of_unity(N, C[0])
roots[C[1]] = util.find_2nth_root_of_unity(N, C[1])

print(roots)

a = RNSPolynomial(B, C, roots, False,[1, 2, 3, 4, 1, 2, 3, 4])#, q=1048721)
b = RNSPolynomial(B, C, roots, False,[1, 0, 0, 0, 0, 0, 0, 0])#, q=1048721)
print(a)
print(b)

#a.convert_RNS_to_NTT()
#b.convert_RNS_to_NTT()

print("a:", a)
print("b:", b)

c = a + b

print("c:", c, "(before rescaling)")
#c.rescale()
c.convert_NTT_to_RNS()
print("c:", c)


exit()
"""
print("setup")
rnsckks = ckks.RNSCKKS(N, p, k, q0, q, L, 2)
print("keygen")
pk, sk = rnsckks.keygen()


m1 = [0.1]*(N//2)
m2 = [0.03]*(N//2)

p1 = rnsckks.encode(m1)
p2 = rnsckks.encode(m2)

"""
c1 = rnsckks.encrypt(p1, pk)
c2 = rnsckks.encrypt(p2, pk)
print("c1:", c1)

c3 = c1 + c2

p3p = rnsckks.decrypt(c3, sk)
m3p = rnsckks.decode(p3p)
print("p3p:", p3p)
print("m3p:",m3p)
"""

print("\n")

print("p1", p1)
print("p2", p2)

p3 = p1 + p1
print("p3", p3)

print("m'1:", rnsckks.decode(p1))
print("m'2:", rnsckks.decode(p2))
print("\n")
print("p3=p1+p1:", rnsckks.decode(p3))

p4 = p3 * p2
p4.rescale()
print("p4=p3*p2:", rnsckks.decode(p4))

p5 = p3 * p1
p5.rescale()
print("p5=p4*p1:", rnsckks.decode(p5))
"""

start = timer()
ckks = ckks.CKKS(N, P, q0, delta, L, 2)
print("Setup duration:", timer() - start)

start = timer()
pk, sk = ckks.keygen()
evk = ckks.evkeygen(sk)
print("Keygen duration:", timer() - start)

#m1 = [0.1, 0.2, 0.3, 0.4]
#m2 = [0.03, 0.04, 1, 1]
m1 = [0.1]*(N//2)
m2 = [0.03]*(N//2)

print("m1: ", m1)
print("m2: ", m2)
print("\n")

start = timer()
p1 = ckks.encode(m1)
print("Encoding duration", timer() - start)
p2 = ckks.encode(m2)
print("p1: ", p1)
print("p2: ", p2)

start = timer()
c1 = ckks.encrypt(p1, pk)
print("Encryption duration", timer() - start)
c2 = ckks.encrypt(p2, pk)

print("c1: ", c1)
print("c2: ", c2)

c_add = c1 + c2
c_add_const = c1 + p2
c_sub = c1 - c2
c_sub_const = c1 - p2
c_mult_const = c1 * p2 * p1
c_mult = c1 * [c_mult_const, evk]

print(c_mult)

print("")
print("e_add", ckks.decode(ckks.decrypt(c_add, sk)))
print("e_add_const", ckks.decode(ckks.decrypt(c_add_const, sk)))
print("e_sub", ckks.decode(ckks.decrypt(c_sub, sk)))
print("e_sub_const", ckks.decode(ckks.decrypt(c_sub_const, sk)))
print("e_mult", ckks.decode(ckks.decrypt(c_mult, sk)))
print("e_mult_const", ckks.decode(ckks.decrypt(c_mult_const, sk)))

"""
