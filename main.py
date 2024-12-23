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
M = 8
delta = 2**8 #128 -> delta
q0 = delta * 2**7
L = 2
#q = #2971215073#8380417#12289#7681#24287#1310717681
#q = 1310717591#8380417

#p = Polynomial([1, 2, 3, 4], 7681)
#print(p)
#print(p+Polynomial([1, 2, 3], 7681))

#pol1 = Polynomial([1, 2, 3, 4], 7681)
#pol2 = Polynomial([5, 6, 7, 8], 7681)

print(util.find_next_prime(2**30))

#print(pol1 * pol2)

ckks = ckks.CKKS(M, util.find_next_prime(q0), util.find_next_prime(delta), L)

m1 = np.array([1+0j, 2+0j])#, 3+0j, 4+0j])
m2 = np.array([3+0j, 4+0j])#, 3+0j, -4+0j])

print("m1: ", m1)
print("m2: ", m2)

print("\n\n")

#p = ckks.encode(np.array([1, 2, 3, 4]))
p1 = ckks.encode(m1)#np.array([1 +5j, 1 - 1j]))
p2 = ckks.encode(m2)#np.array([1 +5j, 2 -1j]))

#comp = 3.+1.j

#print(comp + comp)

#b = [1, 2]#, 3, 4]
#p = ckks.encode(b)
print("p1: ", p1)
print("p2: ", p2)

#h = Polynomial([320, 181, 320, 90], 7681)
#print(ntt.ntt_psi(h))

#print(ntt.ntt_psi(p))

p_add = p1 + p2
print("p1 + p2: ", p_add)

p_mult = p1 * p2
print("p1 * p2: ", p_mult)

p_mult = ckks.rescale(p_mult)
p_mult = p_mult * p_mult
p_mult = ckks.rescale(p_mult)
#p_mult.coeffs = [c / ql for c in p_mult.coeffs]

e_add = ckks.decode(p_add)
e_mult = ckks.decode(p_mult)

#Rescaling test:


#print(p)
print("Result add:  ", e_add)
print("Result mult: ", e_mult)

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