import fourier_transform as fourier
import smwtp
import itertools

instance = smwtp.SMWTP("instances/smtwtp/n5_rdd0.2_tf0.2_seed0.txt")
f = instance.getFunction()
ft = fourier.FourierTransform(5, f, "YKR")

err = 0

for pi in itertools.permutations(range(1,6)):
    print(pi)
    print(smwtp.to_int(pi))
    print(ft.inverseFT(pi, 4))
    print(f[smwtp.to_int(pi)])
    print("--------------------------------")
    err = err + abs(ft.inverseFT(pi, 4) - f[smwtp.to_int(pi)])

err = err / 120
print(err)