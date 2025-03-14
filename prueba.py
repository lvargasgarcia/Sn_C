import fourier_transform
import itertools
from smwtp import *
from functionsamples import *
import random
import scipy
import pandas as pd
import numpy as np

f = {}


instance = FunctionFromSamples("instances/arp/arp_5_17.csv")
dict_instance = instance.getFunction()
# for p in dict_instance:
#     dict_instance[p] = fourier_transform.Rational(int(dict_instance[p]),1)
n = instance.getN()
ft = fourier_transform.FourierTransform_double(n, dict_instance, "YKR", 7)
print("Hola")
ft.set_coefficient(np.ones((5,5)).tolist(), [4,1])

for coef in ft.coefficients.keys():
    print(ft.coefficients[coef])


