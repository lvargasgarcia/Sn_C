import torch
import snob as Snob2
import cnine
from irrep import Irrep
from permutaciones import *
from matrix_utils import *
import unittest
import time
import datetime
import random
import itertools
from fourierTransform import FourierTransform

def log_test_result(test_name, duration):
    with open("resultados.txt", "a") as file:
        file.write(f"{test_name}:{datetime.datetime.now()} passed in {duration:.2f} seconds\n")

dict_gaussian = {tuple(pi): round(random.normalvariate(0, 1), 6) for pi in itertools.permutations([i for i in range(1, 8 + 1)])}



def gaussian(pi):
    return dict_gaussian[tuple(pi)]

class TestIrrepRepresentation(unittest.TestCase):

    def _test_irrep_representation(self, partition, mode="YOR"):
        G = Snob2.Sn(partition.getn())
        rho = Snob2.SnIrrep(partition)
        mi_rho = Irrep(partition, mode)
        if mode != "YOR":
            print("mcm", mi_rho.mcm)

        print("La partición es", partition)

        for i in range(len(G)):
            if i != 0:
                pi = G[i]
                # print("La permutación es", pi)
                # print("La descomposición de la permutación en transposiciones adyacentes es: ", express_into_adyacent_transpositions(pi))

                mi_matriz = mi_rho.evaluate(pi)
                matriz_snob = rho[pi].torch().tolist()

                result = np.allclose(mi_matriz, matriz_snob, atol=1e-6) if mode == "YOR" else True
                self.assertTrue(result, "Las matrices no coinciden")


    def test_YOR_generations(self):

        t_0 = time.time()

        partitions = [[2], [1,1], [3], [2,1], [1,1,1], [4], [3,1], [2,2], [2,1,1], [1,1,1,1],
         [5], [4,1], [3,2], [3,1,1], [2,2,1], [2,1,1,1], [1,1,1,1,1], [6], [5,1],
         [4,2], [4,1,1], [3,2,1], [3,1,1,1], [2,2,1,1], [2,1,1,1,1], [1,1,1,1,1,1],
         [7], [6,1], [5,2], [5,1,1], [4,2,1], [4,1,1,1], [3,2,1,1], [3,1,1,1,1],
         [2,2,1,1,1], [2,1,1,1,1,1], [1,1,1,1,1,1,1], [8], [7,1], [6,2], [6,1,1],
         [5,2,1], [5,1,1,1], [4,2,1,1], [4,1,1,1,1], [3,2,1,1,1], [3,1,1,1,1,1],
         [2,2,1,1,1,1], [2,1,1,1,1,1,1], [1,1,1,1,1,1,1,1], [9], [8,1], [7,2],
         [7,1,1], [6,2,1], [6,1,1,1], [5,2,1,1], [5,1,1,1,1], [4,2,1,1,1],
         [4,1,1,1,1,1], [3,2,1,1,1,1], [3,1,1,1,1,1,1], [2,2,1,1,1,1,1],
         [2,1,1,1,1,1,1,1], [1,1,1,1,1,1,1,1,1], [10], [9,1], [8,2], [8,1,1],
         [7,2,1], [7,1,1,1], [6,2,1,1], [6,1,1,1,1], [5,2,1,1,1], [5,1,1,1,1,1],
         [4,2,1,1,1,1], [4,1,1,1,1,1,1], [3,2,1,1,1,1,1], [3,1,1,1,1,1,1,1],
         [2,2,1,1,1,1,1,1], [2,1,1,1,1,1,1,1,1], [1,1,1,1,1,1,1,1,1,1]]
        
        for partition in partitions:
            partition = Snob2.IntegerPartition(partition)
            n = partition.getn()
            print("n =", )
            print("Probando la partición", partition)
            mi_rho = Irrep(partition, mode="YOR")
            rho = Snob2.SnIrrep(partition)
            for i in range(2,n+1):
                pi = [j for j in range(1,n+1)]
                pi[i-2], pi[i-1] = pi[i-1], pi[i-2]
                print(pi)
                self.assertTrue(np.allclose(mi_rho.evaluate(Snob2.SnElement(pi)), rho[pi].torch().tolist(), atol=1e-6), "Las matrices no coinciden")

        duration = (time.time() - t_0)
        print("Tiempo de ejecución:", duration, "seconds")
        log_test_result("test_YOR_generations", duration)

    def test_YKR_generations(self):

        t_0 = time.time()

        partitions = [[2], [1,1], [3], [2,1], [1,1,1], [4], [3,1], [2,2], [2,1,1], [1,1,1,1],
         [5], [4,1], [3,2], [3,1,1], [2,2,1], [2,1,1,1], [1,1,1,1,1], [6], [5,1],
         [4,2], [4,1,1], [3,2,1], [3,1,1,1], [2,2,1,1], [2,1,1,1,1], [1,1,1,1,1,1],
         [7], [6,1], [5,2], [5,1,1], [4,2,1], [4,1,1,1], [3,2,1,1], [3,1,1,1,1],
         [2,2,1,1,1], [2,1,1,1,1,1], [1,1,1,1,1,1,1], [8], [7,1], [6,2], [6,1,1],
         [5,2,1], [5,1,1,1], [4,2,1,1], [4,1,1,1,1], [3,2,1,1,1], [3,1,1,1,1,1],
         [2,2,1,1,1,1], [2,1,1,1,1,1,1], [1,1,1,1,1,1,1,1], [9], [8,1], [7,2],
         [7,1,1], [6,2,1], [6,1,1,1], [5,2,1,1], [5,1,1,1,1], [4,2,1,1,1],
         [4,1,1,1,1,1], [3,2,1,1,1,1], [3,1,1,1,1,1,1], [2,2,1,1,1,1,1],
         [2,1,1,1,1,1,1,1], [1,1,1,1,1,1,1,1,1], [10], [9,1], [8,2], [8,1,1],
         [7,2,1], [7,1,1,1], [6,2,1,1], [6,1,1,1,1], [5,2,1,1,1], [5,1,1,1,1,1],
         [4,2,1,1,1,1], [4,1,1,1,1,1,1], [3,2,1,1,1,1,1], [3,1,1,1,1,1,1,1],
         [2,2,1,1,1,1,1,1], [2,1,1,1,1,1,1,1,1], [1,1,1,1,1,1,1,1,1,1]]
        
        for partition in partitions:
            partition = Snob2.IntegerPartition(partition)
            print("n =", partition.getn())
            print("Probando la partición", partition)
            mi_rho = Irrep(partition, mode="YKR")
            print("mcm", mi_rho.mcm)
            rho = Snob2.SnIrrep(partition)
        duration = (time.time() - t_0)
        print("Tiempo de ejecución:", duration, "seconds")
        log_test_result("test_YKR_generations", duration)


    def test_partitions_representations(self):

        t_0 = time.time()

        partitions = [[2], [1,1], [3], [2,1], [1,1,1], [4], [3,1], [2,2], [2,1,1], [5], [4,1], [3,2], [3,1,1], [6]
                      , [5,1], [4,2], [4,1,1], [7], [6,1], [5,2], [5,1,1], [8], [7,1], [6,2], [6,1,1]]
        partitions = [Snob2.IntegerPartition(partition) for partition in partitions]
        for partition in partitions:
            print("n =", partition.getn())
            print("Probando la partición", partition)
            self._test_irrep_representation(partition)

        duration = (time.time() - t_0)
        print("Tiempo de ejecución:", duration, "seconds")
        log_test_result("test_partitions_representations", duration)

    def test_partitions_representations_fast(self):

        t_0 = time.time()

        partitions = [[2], [1,1], [3], [2,1], [1,1,1], [4], [3,1], [2,2], [2,1,1], [5], [4,1], [3,2], [3,1,1], [6]
                      , [5,1], [4,2], [4,1,1], [7], [6,1], [5,2], [5,1,1]]
        partitions = [Snob2.IntegerPartition(partition) for partition in partitions]
        for partition in partitions:
            print("n =", partition.getn())
            print("Probando la partición", partition)
            self._test_irrep_representation(partition)

        duration = (time.time() - t_0)
        print("Tiempo de ejecución:", duration, "seconds")
        log_test_result("test_partitions_representations_fast", duration)

    def test_YKR(self):

        t_0 = time.time()

        partitions = [[2], [1,1], [3], [2,1], [1,1,1], [4], [3,1], [2,2], [2,1,1], [5], [4,1], [3,2], [3,1,1], [6]
                      , [5,1], [4,2], [4,1,1], [7], [6,1], [5,2], [5,1,1], [8], [7,1], [6,2], [6,1,1]]
        partitions = [Snob2.IntegerPartition(partition) for partition in partitions]
        for partition in partitions:
            print("n =", partition.getn())
            print("Probando la partición", partition)
            self._test_irrep_representation(partition, mode="YKR")

        duration = (time.time() - t_0)
        print("Tiempo de ejecución:", duration, "seconds")
        log_test_result("test_YKR", duration)

    def test_YSR(self):

        t_0 = time.time()

        partitions = [[2], [1,1], [3], [2,1], [1,1,1], [4], [3,1], [2,2], [2,1,1], [5], [4,1], [3,2], [3,1,1], [6]
                      , [5,1], [4,2], [4,1,1], [7], [6,1], [5,2], [5,1,1], [8], [7,1], [6,2], [6,1,1]]
        partitions = [Snob2.IntegerPartition(partition) for partition in partitions]
        for partition in partitions:
            print("n =", partition.getn())
            print("Probando la partición", partition)
            self._test_irrep_representation(partition, mode="YSR")

        duration = (time.time() - t_0)
        print("Tiempo de ejecución:", duration, "seconds")
        log_test_result("test_YSR", duration)

    def test_YKR_fast(self):

        t_0 = time.time()

        partitions = [[2], [1,1], [3], [2,1], [1,1,1], [4], [3,1], [2,2], [2,1,1], [5], [4,1], [3,2], [3,1,1], [6]
                      , [5,1], [4,2], [4,1,1], [7], [6,1], [5,2], [5,1,1]]
        partitions = [Snob2.IntegerPartition(partition) for partition in partitions]
        for partition in partitions:
            print("n =", partition.getn())
            print("Probando la partición", partition)
            self._test_irrep_representation(partition, mode="YKR")

        duration = (time.time() - t_0)
        print("Tiempo de ejecución:", duration, "seconds")
        log_test_result("test_YKR_fast", duration)

    def test_YSR_fast(self):

        t_0 = time.time()

        partitions = [[2], [1,1], [3], [2,1], [1,1,1], [4], [3,1], [2,2], [2,1,1], [5], [4,1], [3,2], [3,1,1], [6]
                      , [5,1], [4,2], [4,1,1], [7], [6,1], [5,2], [5,1,1]]
        partitions = [Snob2.IntegerPartition(partition) for partition in partitions]
        for partition in partitions:
            print("n =", partition.getn())
            print("Probando la partición", partition)
            self._test_irrep_representation(partition, mode="YSR")

        duration = (time.time() - t_0)
        print("Tiempo de ejecución:", duration, "seconds")
        log_test_result("test_YSR_fast", duration)

    # Test para comprobar que YKR e YSR son representaciones de Sn
    def test_isrep_YKR(self):

        partition = Snob2.IntegerPartition([5,1,1,1,1,1])
        mi_rho = Irrep(partition, mode="YKR")
        G = Snob2.Sn(10)
        for i in range(len(G)):
                if i > 1e2:
                    break
                if i != 0:
                    pi = G[i]
                    matrix = np.array(mi_rho.evaluate(pi), dtype=np.float64)
                    self.assertTrue(np.linalg.det(matrix) != 0, "La matriz de la permutación " +  str(pi) + " no es invertible")
                    print("Es invertible")

    def test_ft_creationtime(self):

        n = 8
        t_0 = time.time()
        ft = FourierTransform(n, gaussian, mode="YKR")
        t_f = time.time()
        print("Tiempo:", t_f - t_0, "s")
        # print(ft)

if __name__ == "__main__":
    unittest.main()

