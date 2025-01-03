from enum import Enum
import random

def gaussian_elimination(matrix, vector):
    n = len(matrix)

    # Ensure the vector supports complex numbers
    vector = [complex(v) for v in vector]

    # Forward elimination
    for i in range(n):
        # Pivot selection (partial pivoting)
        max_row = i + max(range(n - i), key=lambda k: abs(matrix[i + k][i]))
        if matrix[max_row][i] == 0:
            raise ValueError("Matrix is singular or nearly singular.")

        # Swap rows in the matrix and the vector
        matrix[i], matrix[max_row] = matrix[max_row], matrix[i]
        vector[i], vector[max_row] = vector[max_row], vector[i]

        # Make the diagonal element 1 and eliminate below
        for j in range(i + 1, n):
            factor = matrix[j][i] / matrix[i][i]
            for k in range(i, n):
                matrix[j][k] -= factor * matrix[i][k]
            vector[j] -= factor * vector[i]

    # Back substitution
    solution = [0] * n
    for i in range(n - 1, -1, -1):
        sum_ax = sum(matrix[i][j] * solution[j] for j in range(i + 1, n))
        solution[i] = (vector[i] - sum_ax) / matrix[i][i]

    return solution


def findMultInv(a, n):
    # Compute Extended Euclidian algorithm
    # a*x = 1 mod n
    # 1 = a*x + n*y = gcd(a, n)
    gcd, x, _ = extGCD(a, n)
    return x % n

def extGCD(a, b):
    if b == 0:
        return a, 1, 0
    else:
        gcd, x1, y1 = extGCD(b, a % b)
        x = y1
        y = x1 - (a // b) * y1
        return gcd, x, y
        
def GCD(a, b):
    gcd, _, _ = extGCD(a, b)
    return gcd


def totient(m):
    """Calculate Euler's Totient function φ(m) manually."""
    count = 0
    for i in range(1, m):
        if GCD(i, m) == 1:
            count += 1
    return count

def mod_exp(base, exp, mod):
        result = 1
        while exp > 0:
            if exp % 2 == 1:  # If exp is odd
                result = (result * base) % mod
            base = (base * base) % mod
            exp //= 2
        return result

def is_prime(n):
    if n <= 1:
        return False
    if n == 2:
        return True
    if n%2 == 0:
        return False
    
    for i in range(3, int(n**0.5)+1, 2):
        if n % i == 0:
            return False
    return True


def find_next_prime(n):
    p = n
    if is_prime(p):
        return p
    # Go to the next odd number after n
    if n % 2 == 0:
        p+=1
    else:
        p+=2
    
    #Now loop through every odd number to check for primality
    while(True):
        if is_prime(p):
            return p
        p+=2




# Uniformly sample coefficients from Z_q
def sample_uniform_coeffs(n, q):
    return [random.randint(round(-(q-1)/2), round((q-1)/2)) for i in range(0, n)]

# Discrete gaussian distribution with variance sigma^2
def sample_gaussian_coeffs(n, sigma=3.2):
    return [abs(round(random.gauss(0, sigma**2))) for i in range(0, n)]

# Uniform ternany distribution
def sample_uniform_ternary_coeffs(n):
    return [random.choice([-1, 0, 1]) for i in range(0, n)]

