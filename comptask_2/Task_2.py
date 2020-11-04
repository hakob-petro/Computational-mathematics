import numpy as np
from math import sqrt
from math import isnan 

class My_exception(Exception):
    def __init__(self, *args):
        if args:
            text = ""
            for i in range(len(args)):
                text += str(args[i])
            self.text = text
        else:
            self.text = None
            
    def __str__(self):
        return f" My_exception is {self.text}"

    
def square(x):
    return x**2
    
    
def vector_evklid_norm(v):
    return sqrt(sum(map(square, v)))


def discrepancy(b, A, x): #невязка
    return b - mat_mul(A, x)
        
    
def mat_trans(A):
    if A.ndim == 1:
        transposed = A
    else:
        transposed = []
        for row in zip(*list(A)):
            transposed.append(row)
        transposed = np.array(transposed)
    return transposed


def mat_mul(M,P):
    if M.ndim == 1:
        M = np.expand_dims(M, axis=0)
    m, n = M.shape
    
    if P.ndim == 1:
        P = np.expand_dims(P, axis=1)
    p, k = P.shape
    
    if n != p:
        raise My_exception('Wrong dimensions!')
    
    else:
        K = np.zeros((m ,k))
        for i in range(m):
            for j in range(k):
                for l in range(n):
                    K[i][j] += M[i][l]*P[l][j]
    if m == 1:
        K = K[0] #or K.flatten('F')
    if k == 1:
        K = K.flatten()
    return K


def mat_up_triangle(A, N):
    if N == 1:
        up = A
    else:
        up = A.copy()
        for i in range(1, N):
            for j in range(i):
                up[i][:] = up[i][:] - ((up[i][j]/up[j][j])*up[j][:])
    return up

        
def mat_det(A, N):
    A_triag = mat_up_triangle(A, N)
    det = 1
    for i in range(N):
        det *= A_triag[i][i]
    return det
        
def mat_inv(A, N):
    M = np.zeros((N, N))
    for i in range(N):
        help_mat = np.delete(A, i, 0).copy()
        for j in range(N):
            elem = mat_det(np.delete(help_mat, j, 1), N-1)
            if isnan(elem):
                M[i][j] = 0
            else:
                M[i][j] = pow(-1, i+j)*elem
    return (1/mat_det(A, N))*mat_trans(M)


def solv_D_U(A, N):
    D = np.zeros((N,N))
    U = A.copy()
    for i in range(N):
        D[i][i] = A[i][i]
        for j in range(i+1):
            U[i][j] = 0
    return D, U

#==================================================================================================================
#==================================================================================================================

with open('matrix_2.txt') as f:
    matrix = list()
    x = list()
    for row in f.readlines():
        matrix.append(list(map(float, row.split())))
N = int(matrix[0][0])
A = np.array(matrix[1:])
b = np.array([1, 2, 3, 4])#np.random.rand(N)
if sum(map(abs,b)) == 0:
    raise My_exception("b is null vector!")
param_omega = 1.5
epsilon = 0.0001

# Algorithm of successive over-relaxation
D, U = solv_D_U(A, N)
L = A - U - D 

# Condition of positive definite matrix
for i in range(N):
    if mat_det(A[:i+1, :i+1], i) < 0:
        raise My_exception("A isn't a positive definite matrix")
        
# Condition of the symmetry matrix
comparison = U == mat_trans(L)
if not comparison.all():
    raise My_exception("A isn't a symmetry matrix")
    
x_seidel = x_sor = np.random.rand(N) # here x equal to x_0
while vector_evklid_norm(discrepancy(b, A, x_sor)) > epsilon:
    x_seidel = mat_mul(mat_inv(L + D, N),(b - mat_mul(U, x_seidel))) # Seidel L*x_(k+1) + D*x_(k+1) + U*x_k = f 
    x_sor = param_omega * x_seidel + (1 - param_omega) * x_sor # Successive over-relaxation
    print(vector_evklid_norm(discrepancy(b, A, x_sor)))

print(f"The solution is x = {x_sor}")
