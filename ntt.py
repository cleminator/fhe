import math
from poly import Polynomial
import util



def ntt(a, q, root):
    a_hat = []
    n = len(a)
    for j in range(0, n):
        aj = 0
        for i in range(0, n):
            aj += root**(i*j) * a[i]
        a_hat.append(aj % q)
    return a_hat

def intt(a_hat, q, root):
    a = []
    n = len(a_hat)
    inv_n = findMultInv(n, q)
    inv_root = findMultInv(root, q)
    for i in range(0, n):
        ai = 0
        for j in range(0, n):
            ai += inv_root**(i*j) * a_hat[j]
        ai *= inv_n
        a.append(ai % q)
    return a

##################

def ntt_psi(a):#, q, root2n):
    a_hat = []
    root2n = a.find_2nth_root_of_unity(a.n)
    for j in range(0, a.n):
        aj = 0
        for i in range(0, a.n):
            aj += root2n**(2*i*j+i) * a.coeffs[i]
        a_hat.append(aj % a.q)
    return Polynomial(a_hat, a.q)

def intt_psi(a_hat):#, q, root2n):
    #print("START")
    a = []
    root2n = a_hat.find_2nth_root_of_unity(a_hat.n)
    inv_n = util.findMultInv(a_hat.n, a_hat.q)
    #print(inv_n)
    inv_root2n = util.findMultInv(root2n, a_hat.q)
    #print(inv_root2n)
    for i in range(0, a_hat.n):
        ai = 0
        for j in range(0, a_hat.n):
            tmp = util.mod_exp(inv_root2n, (2*i*j+i), a_hat.q)#inv_root2n**(2*i*j+j) % a_hat.q
            ai += (tmp * a_hat.coeffs[j]) % a_hat.q
        ai %= a_hat.q
        ai *= inv_n
        a.append(ai % a_hat.q)
    return Polynomial(a, a_hat.q)

