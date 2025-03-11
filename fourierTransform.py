from queue import PriorityQueue
import json
import numpy as np
from irrep import Irrep
from concurrent.futures import ProcessPoolExecutor
import time
import math
from fractions import Fraction
from multiprocessing import Pool, cpu_count
from gmpy2 import mpq, mpz
from williams_interface import build_coefficient
from permutations_iterator import *

def print_mpq_matrix(matrix):
    return [[str(elem.numerator) + "/" + str(elem.denominator) for elem in row] for row in matrix]

def matmul(A, B):

    # Paso 1: obtener el minimo comun multiplo de los denominadores
    dens_A = [mpz(elem.denominator) for row in A for elem in row]
    dens_B = [mpz(elem.denominator) for row in B for elem in row]

    lcm = np.lcm.reduce(dens_A + dens_B)

    normalized_A = np.array([[mpz(elem) * mpz(lcm // elem.denominator) for elem in row] for row in A], dtype=np.int64)
    normalized_B = np.array([[mpz(elem) * mpz(lcm // elem.denominator) for elem in row] for row in B], dtype=np.int64)

    AB = normalized_A @ normalized_B

    lcm_sqr = mpz(lcm)*mpz(lcm)

    return np.array([[mpq(str(elem) + "/" + str(lcm_sqr)) for elem in row] for row in AB], dtype=object)


def generate_partitions(n):
    if n == 1:
        return [[1]]
    else:
        resp = [[n]]
        for i in range(n-1,0, -1):
            for part in generate_partitions(n-i):
                resp.append([i]+part)
        resp = [sorted(part, reverse=True) for part in resp]
        appeared = set()
        i = 0
        while i < len(resp):
            if tuple(resp[i]) in appeared:
                resp.pop(i)
            else:
                appeared.add(tuple(resp[i]))
                i += 1
        return resp

class FourierTransform:
    
    def __init__(self, n , f, mode="YKR"):

        t_0 = time.time()
        self.n = n
        self.nfact = math.factorial(n)
        #self.f = f
        self.mode = mode
        partitions = generate_partitions(self.n)
        self.images = {}
        inverseFT = [np.array([mpq(0) for _ in range(self.nfact)], dtype=object) for _ in range(n)]

        irreps = {tuple(partition): Irrep(partition, mode=self.mode) for partition in partitions}
        d_lambda = {tuple(partition): irreps[tuple(partition)].matrices[0].shape[0] for partition in partitions}

        for partition in partitions:
            print("d_lambda", partition, "->", d_lambda[tuple(partition)])

        (self.f_nums, self.f_dens, self.williams_sequence) = generate_williams_list(f, n)

        # Crear una PriorityQueue
        priority_queue = PriorityQueue()

        # Insertar las particiones con su prioridad en orden de d_lambda creciente
        for partition in partitions:
            # Prioridad es inversamente proporcional a d_lambda
            priority = d_lambda[tuple(partition)]
            priority_queue.put((priority, partition))
               
        tau = [2,1] + [i for i in range(3,n+1)] # transposición (1,2)
        sigma = [(i+1) for i in range(1,n)] + [1] # Ciclo completo
        q = [(n-i+1) for i in range(1,n+1)] # El 1 va al n, el 2 al n-1, etc.
        qtau = compose(q, tau)
        qsigma = compose(q, sigma)
        qsigmatau = compose(qsigma, tau)
        invSigma = inverse(sigma)

        cpu_cores = cpu_count()
        print("CPU cores:", cpu_cores)

        workers = cpu_cores 
        
        # # Usar ProcessPoolExecutor con la cola de prioridades
        with ProcessPoolExecutor(max_workers=2) as executor:
            
            futures = []
            
            while not priority_queue.empty():
                _, partition = priority_queue.get()
                doInv = 1 
                irrep = irreps[tuple(partition)]
                
                # tau_matrix = irrep.evaluate_floatingpoint(tau)
                # qsigmatau_matrix = irrep.evaluate_floatingpoint(qsigmatau)
                # invSigma_matrix = irrep.evaluate_floatingpoint(invSigma)

                tau_matrix = irrep.evaluate(tau)
                qsigmatau_matrix = irrep.evaluate(qsigmatau)
                invSigma_matrix = irrep.evaluate(invSigma)

                future = executor.submit(self._buildFT, partition, tau_matrix, qsigmatau_matrix, invSigma_matrix, doInv)
                futures.append(future)
            
            for future in futures:
                partition, fourier_coef = future.result()
                ord = self.n - partition[0]
                fourier_coeff, inv_f = fourier_coef
                self.images[partition] = fourier_coeff
                for i in range(n):
                    if i >= ord:
                        inverseFT[i] += inv_f


        # VERSIÓN SIN PARALELIZACIÓN (PARA DEBUGGING)

        # futures = []
        
        # while not priority_queue.empty():
        #     _, partition = priority_queue.get()
        #     doInv = 1 
        #     irrep = irreps[tuple(partition)]
        #     tau_matrix = np.array([[mpq(elem) for elem in row] for row in irrep.evaluate(Snob2.SnElement(tau))], dtype=object)
        #     qsigmatau_matrix = np.array([[mpq(elem) for elem in row] for row in irrep.evaluate(Snob2.SnElement(qsigmatau))], dtype=object)
        #     invSigma_matrix = np.array([[mpq(elem) for elem in row] for row in irrep.evaluate(Snob2.SnElement(invSigma))], dtype=object)
        #     future = self._buildFT(partition, tau_matrix, qsigmatau_matrix, invSigma_matrix, doInv)
        #     futures.append(future)
        
        # for future in futures:
        #     partition, fourier_coef = future
        #     ord = self.n - partition[0]
        #     fourier_coeff, inv_f = fourier_coef
        #     self.images[partition] = fourier_coeff
        #     for i in range(n):
        #         if i >= ord:
        #             inverseFT[i] += inv_f
    
        self.irreps = irreps
        self.d_lambda = d_lambda

        self.inverseFT = {}

        p = compose(qsigma, tau)

        k = 0

        for j in range(n):
            self.inverseFT[(tuple(inverse(p)),(j))] = inverseFT[j][k]

        k += 1

        while(p != qtau):
            
            if(p != qsigmatau):
                if(williamsCondition(p,n) and p != qsigma):
                    p = compose(p, tau)

                else:
                    p = compose(p, invSigma)

            else:
                p = compose(p, invSigma)

            for j in range(n):
                self.inverseFT[(tuple(inverse(p)),(j))] = inverseFT[j][k]

            k += 1

        print("FT construida en", time.time() - t_0, "s")


    def __str__(self):
        np.set_printoptions(threshold=np.inf)
        return json.dumps({str(key): print_mpq_matrix(value) for key, value in self.images.items()})
    
    def _buildFT_williams_YOR(self, tau_matrix, qsigmatau_matrix, invSigma_matrix):

        n = self.n
        tau = [2,1] + [i for i in range(3,n+1)] # transposición (1,2)
        sigma = [(i+1) for i in range(1,n)] + [1] # Ciclo completo
        q = [(n-i+1) for i in range(1,n+1)] # El 1 va al n, el 2 al n-1, etc.
        qsigma = compose(q, sigma)
        qtau = compose(q, tau)
        qsigmatau = compose(qsigma, tau)
        invSigma = inverse(sigma)

        p = compose(qsigma, tau)

        p_matrix = qsigmatau_matrix
        fourier_matrix = self.f(tuple(p))*p_matrix

        while(p != qtau):

            if(p != qsigmatau):
                if(williamsCondition(p,n) and p != qsigma):
                    p = compose(p, tau)
                    p_matrix = p_matrix @ tau_matrix
                    fourier_matrix += self.f(tuple(p))*p_matrix
                else:
                    p = compose(p, invSigma)
                    p_matrix = p_matrix @ invSigma_matrix
                    fourier_matrix += self.f(tuple(p))*p_matrix
            else:
                p = compose(p, invSigma)
                p_matrix = p_matrix @ invSigma_matrix
                fourier_matrix += self.f(tuple(p))*p_matrix
        
        return fourier_matrix

    def _buildFT_williams(self, tau_matrix, qsigmatau_matrix, invSigma_matrix, doInv):

 
    
        # fourier_matrix = mpq(self.f(tuple(p))) * p_matrix

        return build_coefficient(self.n, tau_matrix.shape[0], tau_matrix, qsigmatau_matrix, invSigma_matrix, self.f_nums, self.f_dens, self.williams_sequence, doInv)

        p_matrix = qsigmatau_matrix
        tau_matrix_n = tau_matrix
        invSigma_matrix_n = invSigma_matrix
        i = 0
        fourier_matrix = (self.f_nums[i]/self.f_dens[i]) * p_matrix
        i += 1

        for step in self.williams_sequence:
            if step == 't':
                p_matrix = p_matrix @ tau_matrix_n
            else:
                p_matrix = p_matrix @ invSigma_matrix_n
            fourier_matrix += (self.f_nums[i]/self.f_dens[i]) * p_matrix
            i += 1
            if(tau_matrix.shape[0] > 100):
                print(tau_matrix.shape[0], "->", i)
        
        inverseFT = np.zeros((self.nfact), dtype=np.float64)

        p_matrix = qsigmatau_matrix
        tau_matrix_n = tau_matrix
        invSigma_matrix_n = invSigma_matrix

        inverseFT[0] = np.trace(fourier_matrix @ p_matrix)
        i = 1

        for step in self.williams_sequence:
            if step == 't':
                p_matrix = p_matrix @ tau_matrix_n
            else:
                p_matrix = p_matrix @ invSigma_matrix_n
            inverseFT[i] = np.trace(fourier_matrix @ p_matrix)
            i += 1
            if(tau_matrix.shape[0] > 100):
                print(tau_matrix.shape[0], "->", i)
        
        return fourier_matrix, (tau_matrix.shape[0]/self.nfact)*inverseFT
        
        # while(p != qtau):

        #     if(p != qsigmatau):
        #         if(williamsCondition(p,n) and p != qsigma):
        #             p = compose(p, tau)
        #             # p_matrix = multiplicar_matrices_fraccionadas(p_matrix, tau_matrix)
        #             # p_matrix = multiply_mpq_matrices(p_matrix, tau_matrix)
        #             p_matrix = p_matrix @ tau_matrix

        #             f_p = mpq(self.f(tuple(p)))
        #             fourier_matrix += f_p*p_matrix

        #         else:
        #             p = compose(p, invSigma)
        #             # p_matrix = multiplicar_matrices_fraccionadas(p_matrix, invSigma_matrix)
        #             # p_matrix = multiply_mpq_matrices(p_matrix, invSigma_matrix)
        #             p_matrix = p_matrix @ invSigma_matrix

        #             f_p = mpq(self.f(tuple(p)))
        #             fourier_matrix += f_p*p_matrix

        #     else:
        #         p = compose(p, invSigma)
        #         # p_matrix = multiplicar_matrices_fraccionadas(p_matrix, invSigma_matrix)
        #         # p_matrix = multiply_mpq_matrices(p_matrix, invSigma_matrix)
        #         p_matrix = p_matrix @ invSigma_matrix

        #         f_p = mpq(self.f(tuple(p)))
        #         fourier_matrix += f_p*p_matrix


    
    def _buildFT(self, partition, tau_matrix, qsigmatau_matrix, invSigma_matrix, doInv):
        
        t_0 = time.time()
        if self.mode != "YOR":
            matrix = self._buildFT_williams(tau_matrix, qsigmatau_matrix, invSigma_matrix, doInv)
        else:
            matrix = self._buildFT_williams_YOR(tau_matrix, qsigmatau_matrix, invSigma_matrix)

        print("Construida particion: ", partition, "en", time.time() - t_0, "s")
        
        return tuple(partition), matrix
    
    def _evaluate(self, partition):
        return self.images[tuple(partition)]
    
    def inverseFourierTransform(self, pi, order):
        return self.inverseFT[(pi, (order))]

    # def inverseFourierTransform(self, pi):
    #     f = 0
    #     for partition in self.irreps.keys():
    #         irrep = self.irreps[partition]
    #         matIrrep = irrep.evaluate(Snob2.SnElement(pi).inv())
    #         coefFourier = self.images[partition]
    #         d_lambda = mpq(irrep.matrices[0].shape[0])
    #         f += d_lambda * (np.trace(coefFourier @ matIrrep))
    #     return f/self.nfact
