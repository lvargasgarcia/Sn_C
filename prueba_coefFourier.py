import itertools
import struct
import subprocess
import random
from time import sleep
import torch
import snob as Snob2
import cnine
import sys
import numpy as np
from irrep import Irrep
from fourierTransform import *

def contar_decimales(num):
    # Convertimos el número a cadena, eliminamos el signo "-" si existe
    str_num = str(abs(num))
    
    # Si tiene punto decimal, contamos los dígitos después del punto
    if '.' in str_num:
        return len(str_num.split('.')[1])
    else:
        return 0

def representar_decimal(num):
    decims = contar_decimales(num)
    if decims > 0:
        numerator = int(num * 10**decims)
        denominator = int(10**decims)
        return Fraction(numerator, denominator)
    else:
        return Fraction(num)

def send_info(matrices, im_f_nums, im_f_dens, n):

    process = subprocess.Popen(
        ['valgrind', '--leak-check=full', '--show-leak-kinds=all', '--track-origins=yes', 
'/home/lalo/Sn_representations/williams_walk'],  # El programa C con Valgrind
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE
    )


    process.stdin.write(struct.pack('i', matrices[0].shape[0]))
    for matrix in matrices:
        send_matrix(matrix, process)
    
    im_f_nums = np.array(im_f_nums, dtype=np.int64)
    im_f_dens = np.array(im_f_dens, dtype=np.int64)

    print("python")
    for k in range(len(im_f_nums)):
        print(f"f[{k}]={im_f_nums[k]}/{im_f_dens[k]}")


    process.stdin.write(struct.pack('i', len(im_f_nums)))
    process.stdin.write(struct.pack(f'{im_f_nums.size}q', *im_f_nums))
    process.stdin.write(struct.pack(f'{im_f_dens.size}q', *im_f_dens))
    process.stdin.write(struct.pack('i', n))

    process.stdin.close()

    process.wait()

    print(process.stdout.read().decode('utf-8'))
    print(process.stderr.read().decode('utf-8'))

def send_matrix(matrix, process):
    
    matrix_nums = np.array([mpq(elem).numerator for elem in matrix.flatten()], dtype=np.int64)
    matrix_dens = np.array([mpq(elem).denominator for elem in matrix.flatten()], dtype=np.int64)

    # Enviar la matriz (todos los valores como tipo 'double')
    matrix_bytes = struct.pack(f'{matrix_nums.size}q', *matrix_nums.flatten())
    process.stdin.write(matrix_bytes)

    matrix_bytes = struct.pack(f'{matrix_dens.size}q', *matrix_dens.flatten())
    process.stdin.write(matrix_bytes)


partition = [2,2,1]

n = sum(partition)

irrep = Irrep(Snob2.IntegerPartition(partition), mode="YKR")

tau = [2,1] + [i for i in range(3,n+1)] # transposición (1,2)
sigma = [(i+1) for i in range(1,n)] + [1] # Ciclo completo
q = [(n-i+1) for i in range(1,n+1)] # El 1 va al n, el 2 al n-1, etc.
qtau = compose(q, tau)
qsigma = compose(q, sigma)
qsigmatau = compose(qsigma, tau)
invSigma = inverse(sigma)

tau_matrix = np.array([[mpq(elem) for elem in row] for row in irrep.evaluate(Snob2.SnElement(tau))], dtype=object)
qsigmatau_matrix = np.array([[mpq(elem) for elem in row] for row in irrep.evaluate(Snob2.SnElement(qsigmatau))], dtype=object)
invSigma_matrix = np.array([[mpq(elem) for elem in row] for row in irrep.evaluate(Snob2.SnElement(invSigma))], dtype=object)
p_matrix = qsigmatau_matrix

dim = p_matrix.shape[0]

matrices = [tau_matrix, qsigmatau_matrix, invSigma_matrix]

dict_gaussian = {tuple(p): representar_decimal(round(random.normalvariate(0,1), 2)) for p in itertools.permutations(range(1, n + 1))}

def f(p):
    return dict_gaussian[p]

im_f_nums = []
im_f_dens = []

p = compose(qsigma, tau)
val = mpq(f(tuple(p)))
im_f_nums.append(val.numerator)
im_f_dens.append(val.denominator)

while(p != qtau):

    if(p != qsigmatau):
        if(williamsCondition(p,n) and p != qsigma):
            p = compose(p, tau)
            val = mpq(f(tuple(p)))
            im_f_nums.append(val.numerator)
            im_f_dens.append(val.denominator)
        else:
            p = compose(p, invSigma)
            val = mpq(f(tuple(p)))
            im_f_nums.append(val.numerator)
            im_f_dens.append(val.denominator)
    else:
        p = compose(p, invSigma)
        val = mpq(f(tuple(p)))
        im_f_nums.append(val.numerator)
        im_f_dens.append(val.denominator)

send_info(matrices, im_f_nums, im_f_dens, n)

# p = compose(qsigma, tau)
