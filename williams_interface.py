from permutations_iterator import generate_williams_list
import numpy as np
import math
import ctypes
from fractions import Fraction
from gmpy2 import mpq

class MPQ(ctypes.Structure):
    _fields_ = [("num", ctypes.c_void_p),  # Pointer to numerator (mpz_t)
                ("den", ctypes.c_void_p)]  # Pointer to denominator (mpz_t)]

def build_coefficient(n, dim, tau_matrix, qsigmatau_matrix, invSigma_matrix, f_nums, f_dens, williams_sequence, doInv):
    # Cargar la biblioteca compilada (.so)
    lib = ctypes.CDLL("./libwilliams.so")

    # Definir la función en ctypes
    lib.williams_wrapper.argtypes = [
        ctypes.c_int, ctypes.c_int, ctypes.c_int,  # n, nfact, dim
        ctypes.POINTER(ctypes.c_long), ctypes.POINTER(ctypes.c_long),  # tau_nums, tau_dens
        ctypes.POINTER(ctypes.c_long), ctypes.POINTER(ctypes.c_long),  # qsigmatau_nums, qsigmatau_dens
        ctypes.POINTER(ctypes.c_long), ctypes.POINTER(ctypes.c_long),  # invSigma_nums, invSigma_dens
        ctypes.POINTER(ctypes.c_long), ctypes.POINTER(ctypes.c_long),
        ctypes.c_char_p, ctypes.c_int  # williams_sequence
        # ctypes.POINTER(ctypes.c_long), ctypes.POINTER(ctypes.c_long)  # f_nums, f_dens
    ]
    lib.williams_wrapper.restype = ctypes.c_char_p 

    # Definir constantes
    nfact = math.factorial(n)

    tau_nums = np.array([[elem.numerator for elem in row] for row in tau_matrix], dtype=np.int64)
    tau_dens = np.array([[elem.denominator for elem in row] for row in tau_matrix], dtype=np.int64)
    qsigmatau_nums = np.array([[elem.numerator for elem in row] for row in qsigmatau_matrix], dtype=np.int64)
    qsigmatau_dens = np.array([[elem.denominator for elem in row] for row in qsigmatau_matrix], dtype=np.int64)
    invSigma_nums = np.array([[elem.numerator for elem in row] for row in invSigma_matrix], dtype=np.int64)
    invSigma_dens = np.array([[elem.denominator for elem in row] for row in invSigma_matrix], dtype=np.int64)

    resp = lib.williams_wrapper(
        n, nfact, dim,
        tau_nums.ctypes.data_as(ctypes.POINTER(ctypes.c_long)),
        tau_dens.ctypes.data_as(ctypes.POINTER(ctypes.c_long)),
        qsigmatau_nums.ctypes.data_as(ctypes.POINTER(ctypes.c_long)),
        qsigmatau_dens.ctypes.data_as(ctypes.POINTER(ctypes.c_long)),
        invSigma_nums.ctypes.data_as(ctypes.POINTER(ctypes.c_long)),
        invSigma_dens.ctypes.data_as(ctypes.POINTER(ctypes.c_long)),
        f_nums.ctypes.data_as(ctypes.POINTER(ctypes.c_long)),
        f_dens.ctypes.data_as(ctypes.POINTER(ctypes.c_long)),
        williams_sequence.encode('utf-8'),
        doInv
    )

    if doInv == 0:
    
        resp = resp.decode('utf-8').split(" ")
        aux = np.empty((dim, dim), dtype=object)
        for i in range(dim):
            for j in range(dim):
                aux[i][j] = mpq(resp[i*dim + j])
        
        return aux
    
    else:
    
        resp = resp.decode('utf-8').split("$")
        resp1 = resp[0].split(" ")
        resp2 = resp[1].split(" ")
        aux1 = np.empty((dim, dim), dtype=object)
        aux2 = np.empty((nfact), dtype=object)
        for i in range(dim):
            for j in range(dim):
                aux1[i][j] = mpq(resp1[i*dim + j])
        for i in range(nfact):
            aux2[i] = mpq(resp2[i])
        
        return aux1, aux2




