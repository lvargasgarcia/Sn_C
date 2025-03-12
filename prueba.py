import fourier_transform
import itertools
from smwtp import *
from functionsamples import *
import random 

f = {}

instance = FunctionFromSamples("instances/arp/arp_5_17.csv")
dict_instance = instance.getFunction()
for p in dict_instance:
    dict_instance[p] = fourier_transform.Rational(int(dict_instance[p]),1)
n = instance.getN()
ft = fourier_transform.FourierTransform_Rational(n, dict_instance, "YKR", 6)


for pi in itertools.permutations([i for i in range(1,n+1)]):
    print(dict_instance[to_int(pi)])
    print(ft.inverseFT(pi,n-1))
    print("-------------------")