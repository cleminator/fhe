import math
#from poly import Polynomial
import util

# Schoolbook version implemented after https://eprint.iacr.org/2024/585.pdf

"""
def ntt(a):
    a_hat = []
    root = a.find_primitive_nth_root_of_unity(a.n)
    for j in range(0, a.n):
        aj = 0
        for i in range(0, a.n):
            aj += root**(i*j) * a.coeffs[i]
        a_hat.append(aj % a.q)
    return Polynomial(a_hat, a.q)

def intt(a_hat):
    a = []
    root = a_hat.find_primitive_nth_root_of_unity(a_hat.n)
    inv_n = util.findMultInv(a_hat.n, a_hat.q)
    inv_root = util.findMultInv(root, a_hat.q)
    for i in range(0, a_hat.n):
        ai = 0
        for j in range(0, a_hat.n):
            ai += inv_root**(i*j) * a_hat.coeffs[j]
        ai *= inv_n
        a.append(ai % a_hat.q)
    return Polynomial(a, a_hat.q)
"""
##################

def ntt_psi(coeffs, n, q, root2n):#, q, root2n):
    a_hat = []

    for j in range(0, n):
        #print("j:", j)
        aj = 0
        for i in range(0, n):
            #print("i:", i)
            aj += util.mod_exp(root2n, 2*i*j+i, q) * coeffs[i]#root2n**(2*i*j+i) * coeffs[i]
        a_hat.append(util.mod(aj, q))
    return a_hat#Polynomial(a_hat, q)

def intt_psi(coeffs, n, q, root2n):#, q, root2n):
    #print("START")
    a = []
    inv_n = util.findMultInv(n, q)
    inv_root2n = util.findMultInv(root2n, q)

    for i in range(0, n):
        #print("i:", i)
        ai = 0
        for j in range(0, n):
            #print("j:", j)
            tmp = util.mod_exp(inv_root2n, 2*i*j+i, q)
            ai += util.mod((tmp * coeffs[j]), q)
        ai %= q
        ai *= inv_n
        a.append(ai % q)
    return a#Polynomial(a, a_hat.q)

