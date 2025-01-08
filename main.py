from poly import Polynomial, RNSPolynomial
import ckks
import util


N = 2**2
delta = 2**40
q0 = delta * 2**7
L = 3
P = 2**90

B = [5]
#C = [13, 17, 19, 23]
C = [1313041, 131071, 131063]
#131072

#131041
#131059
#131063
#131071
#131101
#131111

rnsckks = ckks.RNSCKKS(N, B, C, 2**17, 2)
m1 = [0.1]*(N//2)
m2 = [0.03]*(N//2)

p1 = rnsckks.encode(m1)
p2 = rnsckks.encode(m2)
#p1 = RNSPolynomial(B, C, [-5, 15, 17, -1])
#p2 = RNSPolynomial(B, C, [10, 20, 30, 40])

print(p1)
print(p2)

p3 = p1*p2

print("p1:", rnsckks.decode(p1))
print("p2:", rnsckks.decode(p2))
p3.rescale()
print("p3=p1*p2:", rnsckks.decode(p3))

p4 = p3 * p1
p4.rescale()
print("p4=p3*p1:", rnsckks.decode(p4))

exit()

ckks = ckks.CKKS(N, P, q0, delta, L, 2)

pk, sk = ckks.keygen()
evk = ckks.evkeygen(sk)

#m1 = [0.1, 0.2, 0.3, 0.4]
#m2 = [0.03, 0.04, 1, 1]
m1 = [0.1]*(N//2)
m2 = [0.03]*(N//2)

print("m1: ", m1)
print("m2: ", m2)
print("\n")

p1 = ckks.encode(m1)
p2 = ckks.encode(m2)
print("p1: ", p1)
print("p2: ", p2)

c1 = ckks.encrypt(p1, pk)
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

