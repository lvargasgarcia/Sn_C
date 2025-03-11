import torch
import snob as Snob2
import cnine
import sys
import numpy as np
import math
import argparse
from fractions import Fraction
from irrep import Irrep
from snob_base import SnFunction
import itertools
import random
import time
from fourierTransform import FourierTransform
from permutations_iterator import representar_decimal 
from gmpy2 import mpq
import matplotlib.pyplot as plt

# Configurar argparse para recibir n como argumento
def parse_args():
    parser = argparse.ArgumentParser(description="Ejecutar Fourier Transform con truncamiento.")
    parser.add_argument("n", type=int, help="Número de elementos en la permutación")
    return parser.parse_args()

# Obtener argumentos desde línea de comandos
args = parse_args()
n = args.n

random.seed(44)

dict_gaussian = {pi: round(random.uniform(10,50), 10) for pi in itertools.permutations([i for i in range(1, n + 1)])}

def gaussian(pi):
    return dict_gaussian[tuple(pi)]

mode="YSR"

# ks = [i for i in range(n)]
# max_errores = []

# t_inicial = time.time()

# for k in ks:
#     t_0 = time.time()
#     ft = FourierTransform(n, gaussian, mode=mode, truncating_order=k)
#     t_f = time.time()
#     tiempo_ft = t_f - t_0

#     print("Tiempo truncando en orden:", k, tiempo_ft)

#     errores = []
#     for pi in itertools.permutations([i for i in range(1, n + 1)]):
#         gauss = gaussian(pi)
#         gauss = mpq(representar_decimal(gauss))
#         inv_ft = ft.inverseFT[pi]
#         errores.append(abs(gauss - inv_ft))
    
#     errormax = max(errores).numerator / max(errores).denominator
#     max_errores.append(errormax)
#     print("Error máximo, orden:", k, errormax)

# # # Mejorar el gráfico
# # plt.figure(figsize=(10, 6))
# # plt.plot(ks, max_errores, marker='o', linestyle='-', color='b', label='Error Máximo')
# # plt.xlabel('Orden de Truncamiento')
# # plt.ylabel('Error Máximo')
# # plt.title('Error Máximo vs Orden de Truncamiento')
# # plt.legend()
# # plt.grid(True)
# # plt.savefig('error_maximo_vs_orden_truncamiento.png')

# t_final = time.time()
# tiempo_total = t_final - t_inicial
# print("Tiempo total:", tiempo_total)

ft = FourierTransform(n, gaussian, mode=mode, truncating_order=None)

print(ft)