def gaussian_elimination(matrix, vector):
    
    n = len(matrix)

    # Forward elimination
    for i in range(n):
        # Make the diagonal element non-zero by swapping rows if necessary
        if matrix[i][i] == 0:
            for j in range(i + 1, n):
                if matrix[j][i] != 0:
                    matrix[i], matrix[j] = matrix[j], matrix[i]
                    vector[i], vector[j] = vector[j], vector[i]
                    break
            else:
                raise ValueError("Matrix is singular or nearly singular.")
        
        # Scale the pivot row
        pivot = matrix[i][i]
        for k in range(i, n):
            matrix[i][k] /= pivot
        vector[i] /= pivot

        # Eliminate the current column in rows below
        for j in range(i + 1, n):
            factor = matrix[j][i]
            for k in range(i, n):
                matrix[j][k] -= factor * matrix[i][k]
            vector[j] -= factor * vector[i]

    # Back substitution
    solution = [0] * n
    for i in range(n - 1, -1, -1):
        solution[i] = vector[i]
        for j in range(i + 1, n):
            solution[i] -= matrix[i][j] * solution[j]

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

