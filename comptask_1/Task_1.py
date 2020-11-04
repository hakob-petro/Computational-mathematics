#Computer task №1
#Petrosyan Hakob, group 814
"""
Disclaimer!
In the code, I use matrices to further use the full power of numpy.

P.S. I have hardly used any functions from numpy, except for the inverse matrix
"""

import numpy as np
from numpy import linalg as LA
from math import sqrt

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
        print("Wrong dimensions!")
        raise SystemExit
    
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

#==================================================================================================
#Solving L
with open('matrix_spd.txt') as f:
    matrix = list()
    x = list()
    for row in f.readlines():
        matrix.append(list(map(float, row.split())))
N = int(matrix[0][0])
for i in range(N):
    x.append(i+1)
x = np.array(x)
A = np.array(matrix[1:])
b = mat_mul(A, x)
L = np.zeros((N,N))


#first step of finding elemts of matrix L
try:
    L[0][0] = sqrt(A[0][0])
except ValueError:
    print("Извлечение корня из отрицательного числа")
    raise SystemExit    
     
    
#second step of finding elemts of matrix L
try:
    for j in range(N-1):
        L[j+1][0] = A[j+1][0] / L[0][0]
except ZeroDivisionError:
    print("Деление на ноль")
    raise SystemExit


#third step of finding elemts of matrix L
try:
    for i in range(1, N):
        for j in range(i, N):
            summ_ii = 0
            summ_ji = 0
        
            if j == i:
                for p in range(i):
                    summ_ii += L[i][p]**2 
                L[i][i] = sqrt(A[i][i] - summ_ii)
                summ_ii = 0
            
            if j != i:
                for p in range(i):
                    summ_ji += L[i][p]*L[j][p]
                L[j][i] = (A[j][i] - summ_ji) / L[i][i]
                summ_ji = 0
        
except ZeroDivisionError:
    print("Деление на ноль")
    raise SystemExit
except ValueError:
    print("Извлечение корня из отрицательного числа")
    raise SystemExit
#==================================================================================================  
    
#Algorihm
L.round(2)
L_inv = LA.inv(L)
L_T_inv = LA.inv(mat_trans(L))
x_help =  (mat_mul(L_inv, b)).round(2)
x_final = (mat_mul(L_T_inv, x_help)).round(2)
print(x_final - x)
print(x_final)
