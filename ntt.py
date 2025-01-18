import math
import util

# Schoolbook version implemented after https://eprint.iacr.org/2024/585.pdf

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

#############################

def bit_reverse_copy(a):

    n = len(a)
    A = [0] * n  # Initialize an array of size n with zeros

    def rev(k, bits):
        """Compute the bit-reversed value of k for a given number of bits."""
        result = 0
        for _ in range(bits):
            result = (result << 1) | (k & 1)
            k >>= 1
        return result

    bits = n.bit_length() - 1  # Number of bits needed for indices in [0, n-1]
    for k in range(n):
        A[rev(k, bits)] = a[k]

    return A

def fast_ntt(a, n, q, root2n, s=1):
    # Using Cooley-Tukey Algorithm (taken from https://eprint.iacr.org/2024/585.pdf)

    if n == 1:
        return a
    else:
        #a_hat = []
        root_sq = util.mod_exp(root2n, 2, q)
        a_even = a[::2]
        a_odd = a[1::2]

        b_even = fast_ntt(a_even, n//2, q, root_sq, 2*s)
        b_odd = fast_ntt(a_odd, n//2, q, root_sq, 2*s)

        b = [0]*n
        #rk = 1

        for k in range(n // 2):
            rk = util.mod_exp(root2n, 2*k+1, q)
            t = util.mod(rk * b_odd[k], q)
            b[k] = util.mod(b_even[k] + t, q)
            b[k + n//2] = util.mod(b_even[k] - t, q)

            #rk = util.mod(rk * root2n, q)

        return b




