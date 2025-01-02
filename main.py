import ntt
from poly import Polynomial
import ckks
import util
import numpy as np
"""
g = Polynomial([1, 2, 3, 4], 7681)
h = Polynomial([320, 181, 320, 90], 7681)
print("g & h")
g_hat = ntt.ntt_psi(g)
print(g_hat)
h_hat = ntt.ntt_psi(h)
print(h_hat)
print(g * h)

print("\n\n...\n\n")

g_hat = ntt.ntt_psi(g)
print(g_hat)

print(ntt.intt_psi(g_hat))

h_hat = ntt.ntt_psi(h)
print(h_hat)
print(ntt.intt_psi(h_hat))
"""
print(util.sample_gaussian_coeffs(10))


#exit()

M = 8
delta = 2**30 #128 -> delta
q0 = delta * 2**6
L = 2

#q = #2971215073#8380417#12289#7681#24287#1310717681
#q = 1310717591#8380417

#p = Polynomial([1, 2, 3, 4], 7681)
#print(p)
#print(p+Polynomial([1, 2, 3], 7681))

#pol1 = Polynomial([1, 2, 3, 4], 7681)
#pol2 = Polynomial([5, 6, 7, 8], 7681)

#print(util.find_next_prime(2**30))

#print(pol1 * pol2)

ckks = ckks.CKKS(M, util.find_next_prime(q0), util.find_next_prime(delta), L, 2)

pk, sk = ckks.keygen()
evk = ckks.evkeygen(sk)

m1 = np.array([0.1, 0.2])#, 3+0j, 4+0j])
m2 = np.array([0.03, 0.04])#, 3+0j, -4+0j])

print("m1: ", m1)
print("m2: ", m2)


#p = ckks.encode(np.array([1, 2, 3, 4]))
p1 = ckks.encode(m1)#np.array([1 +5j, 1 - 1j]))
p2 = ckks.encode(m2)#np.array([1 +5j, 2 -1j]))

#print(type(p1), isinstance(p1, Polynomial))
#print(type(-1), isinstance(1, Polynomial))

#comp = 3.+1.j

#print(comp + comp)

#b = [1, 2]#, 3, 4]
#p = ckks.encode(b)
print("p1: ", p1)
print("p2: ", p2)

#h = Polynomial([320, 181, 320, 90], 7681)
#print(ntt.ntt_psi(h))

#print(ntt.ntt_psi(p))
#e = p1*p2
#e.rescale(delta)
#print(ckks.decode(e))



c1 = ckks.encrypt(p1, pk)
#print("C1")
#print(p1)
#print(c1.b)
#print(c1.a)
#print(pk[0])
#print(pk[1])

c2 = ckks.encrypt(p2, pk)
print("c1: ", c1)
print("Addition: ")
c3 = c2*[c1, evk]#c1 * p1
print("c3: ", c3)
#c3 = c3 - c2
res = ckks.decrypt(c3, sk)


#print("pp1: ", pp1)
#print("ep1: ", ckks.decode(pp1))
print("e1:  ", ckks.decode(res))

exit()

p_add = p1 + p2
print("p1 + p2: ", p_add)
p_add2 = -1 * p_add
print("scalar: ", p_add2)

p_mult = p1 * (p2, p2)
print("p1 * p2: ", p_mult)

#p_mult.rescale(delta)

#p_mult = p1 * p_mult# * p1
#p_mult = ckks.rescale(p_mult)
#p_mult.coeffs = [c / ql for c in p_mult.coeffs]

print(p_mult)

e_add = ckks.decode(p_add)
e_mult = ckks.decode(p_mult)
e_add2 = ckks.decode(p_add2)

#Rescaling test:


#print(p)
print("Result add:  ", e_add)
print("Result mult: ", e_mult)
print("Result add * 2: ", e_add2)

#print("\n\n")

#t1 = np.polynomial.Polynomial([768, -181, 0, 181][::-1])
#t2 = np.polynomial.Polynomial([-1443, 2500, 0, -2500][::-1])
#mo = np.polynomial.Polynomial([1, 0, 0, 0, 1])

#res1 = t1 * t2 % mo

#print(res1)
#res = p + p

#res2 = Polynomial([0, -2181183, 2013224, -2181183][::-1], 7681)

#print(ckks.decode(res2))



"""
print(ntt.findMultInv(4, 7681))

g = [1, 2, 3, 4]
h = [5, 6, 7, 8]
q = 7681
n = 4
root = ntt.find_nth_root_of_unity(n, q)#3383
root2n = ntt.find_2nth_root_of_unity(n, q)

print(root2n, ntt.findMultInv(root2n, q), ntt.findMultInv(n, q))

#print(root, ", ", root2n)

print(ntt.ntt_psi(h, q, root2n))

print(ntt.intt_psi([2489, 7489, 6478, 6607], q, root2n))

"""