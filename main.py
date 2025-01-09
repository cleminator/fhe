from poly import Polynomial, RNSPolynomial
import ckks
from ntt import ntt, intt, ntt_psi, intt_psi
import util
from timeit import default_timer as timer

N = 2**3
q = 20
q0 = 30
L = 2
#P = 2**90
p = 30


a = Polynomial([1, 2, 3, 4, 1, 2, 3, 4], q=1048721)
b = Polynomial([1, 0, 0, 1, 0, 0, 0, 0], q=1048721)
root2n = a.find_2nth_root_of_unity(a.n)

c = a - b

print(a)
print(b)
print(c)
print("---")
a_hat = ntt_psi(a, root2n)
print(intt_psi(a_hat, root2n))
b_hat = ntt_psi(b, root2n)
print(intt_psi(b_hat, root2n))

mult_d = [ac - bc for ac, bc in zip(a_hat.coeffs, b_hat.coeffs)]
d = Polynomial(mult_d, q=1048721)
print(intt_psi(d, root2n))

exit()
rnsckks = ckks.RNSCKKS(N, p, q0, q, L, 2)
m1 = [0.1]*(N//2)
m2 = [0.03]*(N//2)



p1 = rnsckks.encode(m1)
p2 = rnsckks.encode(m2)

print("p1, p2")
print(p1)
print(p2)

p3 = p1 + p1


print("m'1:", rnsckks.decode(p1))
print("m'2:", rnsckks.decode(p2))
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
