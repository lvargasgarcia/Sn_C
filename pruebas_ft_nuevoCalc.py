import numpy as np
import math
import itertools
import random
import time
from fourierTransform import FourierTransform
from permutations_iterator import representar_decimal 
from gmpy2 import mpq
import matplotlib.pyplot as plt

n = 7

random.seed(44)

dict_gaussian = {pi: (pi[0]) for pi in itertools.permutations([i for i in range(1, n + 1)])}

def gaussian(pi):
    return dict_gaussian[tuple(pi)]

for pi in itertools.permutations([i for i in range(1, n + 1)]):
    print(pi, "-->", gaussian(pi))

mode="YKR"

t_0 = time.time()
ft = FourierTransform(n, gaussian, mode=mode)
t_f = time.time()
fmax = None
fmin = None

tiempo_ft = t_f - t_0
NMAES = []

for k in range(n):
    NMAE = []

    for pi in itertools.permutations([i for i in range(1, n + 1)]):
        gauss = mpq(representar_decimal(gaussian(pi)))

        if fmax is None:
            fmax = gauss
            fmin = gauss
        else:
            if gauss > fmax:
                fmax = gauss
            if gauss < fmin:
                fmin = gauss

        inv_ft = ft.inverseFourierTransform(pi, order=k)
        # print("función", gauss)
        # print("Transformada de Fourier inversa:", inv_ft)
        NMAE.append(abs(gauss - inv_ft))

    NMAE = sum(NMAE)/(math.factorial(n)*abs(fmax - fmin))
    NMAES.append(NMAE)

print(NMAES)

plt.plot(range(n), NMAES)
plt.grid()
plt.xlabel("Orden de truncamiento")
plt.ylabel("NMAE")
#plt.yscale("log")
plt.title("NMAE vs. Orden de truncamiento")
plt.show()