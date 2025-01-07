from poly import Polynomial
import ckks
import util


M = 2**3#32
P = 2**90
delta = 2**30 #128 -> delta
q0 = delta * 2**6
L = 3

ckks = ckks.CKKS(M, P, q0, delta, L, 2)

pk, sk = ckks.keygen()
evk = ckks.evkeygen(sk)

#m1 = [0.1, 0.2, 0.3, 0.4]
#m2 = [0.03, 0.04, 1, 1]
m1 = [0.1]*(M//4)
m2 = [0.03]*(M//4)

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

