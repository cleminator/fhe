from poly import Polynomial, RNSPolynomial
import ckks
from ntt import ntt_psi, intt_psi, fast_ntt, fast_intt
import util
from timeit import default_timer as timer



N = 2**6
q = 30
q0 = 40
L = 2
k = 1
#P = 2**90
p = 18


print("Setup")
rnsckks = ckks.RNSCKKS(N, p, k, q0, q, L, 2)
print("Keygen")
pk, sk = rnsckks.keygen()


m1 = [0.1]*(N//2)
m2 = [0.03]*(N//2)

print("Encoding m1")
p1 = rnsckks.encode(m1)
print("Encoding m2")
p2 = rnsckks.encode(m2)


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
print("p5=p3*p1:", rnsckks.decode(p5))

#c1 = rnsckks.encrypt(p1, pk)
#c2 = rnsckks.encrypt(p2, pk)
#print("\nc1:", c1)
#print("c2:", c2)
#print("m1'", rnsckks.decode(rnsckks.decrypt(c1, sk)))
#print("m2'", rnsckks.decode(rnsckks.decrypt(c2, sk)))

"""

N = 2**3
delta = 2**30
q0 = 2**40
L = 4
P = 2**90

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

