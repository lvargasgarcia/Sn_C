import numpy as np
import math
from gmpy2 import mpq
from permutaciones import express_into_adyacent_transpositions
from matrix_utils import decompress_YOR, decompress, print_matrix
from youngTableau import *

# Aquí se implementará el objeto que representa una representación irreducible de una partición alpha

# El objeto consisirá en un array de n elementos, (alpha es una partición de n)

# Cada elemento será una matriz comprimida que contenga la info de rho^alpha(t_i), estos elementos pueden estar ya calculados o no

# Se irán calculando según se necesiten

# Además, se implementará una función evaluate que dada una permutación, la descomponga en ciclos adyacentes, tome las matrices que necesite y las multiplique

# Eso nos dará una matriz de la representación irreducible de la partición alpha en cualquier elemento de S_n

class Irrep:
    
    # Constructor de la clase, en principio, el array de matrices estará vacío y se irán calculando según haga falta

    def __init__(self, partition, mode: str):
        """
        Inicializa una nueva instancia de IrreducibleRepresentation.

        Args:
            partition (Snob2.IntegerPartition): La partición entera asociada.
            mode (str): El modo de la representación. Debe ser uno de 'YKR', 'YSR' o 'YOR'.

        Raises:
            ValueError: Si el modo no es uno de 'YKR', 'YSR' o 'YOR'.
        """
        if mode not in ["YKR", "YSR", "YOR"]:
            raise ValueError("Mode must be 'YKR', 'YSR', or 'YOR'")
        self.partition = partition
        self.n = sum(partition)
        self.mode = mode
        self.matrices = [None for _ in range(self.n - 1)]
        for i in range(2, self.n + 1):
            self._calc_matrix(i)
    
    def _search_young_tableau(self, tableau):
        shape = tableau.shape()
        n = shape.getn()
        n1 = (-1,-1)
        n2 = (-1,-1)
        
        for i in range(n):
            for j in range(shape[i]):
                if tableau[i, j] == n:
                    n1 = (i,j)
                if tableau[i, j] == n - 1:
                    n2 = (i,j)
                if n1 != (-1,-1) and n2 != (-1,-1):
                    return (n1, n2)  # Devolver tuplas

        return None  # Si no se encuentra, devolver None

    def _n_pos(self, tableau):
        shape = tableau.shape()
        height = shape.height()
        n = shape.getn()
        for i in range(height):
            for j in range(shape[i]):
                if tableau[i, j] == n:
                    n1 = i
                    return n1
        return None 

    def _get_partition_from_gamma(self, y_symbol, partition):
        resp = [] 

        for i in range(len(partition)):
            k = partition[i]
            (a, b) = y_symbol # Desempaquetar la tupla
            if a == i:
                k = k - 1
            if b == i:
                k = k - 1
            resp.append(k)
        while resp and resp[-1] == 0:
            resp.pop()

        return resp

    # ahora vamos a implementar el calculo de la fórmula de frame-robinson-thrall para una particion 

    def _cells_below_cell(self, partition, i, j):
        resp = 0
        for k in range(i + 1, len(partition)):
            if(partition[k] >= j + 1):
                resp = resp + 1
            else:
                break
        return resp

    def _hook_length_prod(self, partition):
        resp = 1
        for i in range(len(partition)):
            for j in range(partition[i]):
                cells_in_row = (partition[i] - j)
                cells_below = self._cells_below_cell(partition, i, j)
                hook_length = (cells_in_row + cells_below)
                resp = resp * hook_length
        return resp

    def _frame_robinson_thrall(self, partition):
        n = sum(partition)
        n_factorial = math.factorial(n)
        den = self._hook_length_prod(partition)
        return int(n_factorial / den)
    
    def _generate_representation_matrix(self, partition, ordered_2ndlev, mode="YKR"):

        n = len(ordered_2ndlev)
        matrix = [[[-1,-1,None] for _ in range(n)] for _ in range(n)]
        
        # Primero hacemos la diagonal
        
        for i in range(n):

            ((a,b),(c,d)) = ordered_2ndlev[i]
            gamma = self._get_partition_from_gamma((a,c), partition)
            d_gamma = self._frame_robinson_thrall(gamma)
            coef = None
            # Cláusula 3 del teorema de young
            if a == c: # R-chain
                coef = 1
            if b == d: # C-chain
                coef = -1
            matrix[i][i] = [d_gamma, d_gamma, coef]

        for i in range(n):
            for j in range(n):
                
                if i < j:

                    ((a,b),(c,d)) = ordered_2ndlev[i]
                    ((a2,_),(c2,_)) = ordered_2ndlev[j]
                    gamma = self._get_partition_from_gamma((a,c), partition)
                    gamma2 = self._get_partition_from_gamma((a2,c2), partition)
                    d_gamma = matrix[i][i][0]
                    d_gamma2 = matrix[j][j][0]
                    
                    matrix[i][j][0] = d_gamma
                    matrix[i][j][1] = d_gamma2

                    matrix[j][i][0] = d_gamma2
                    matrix[j][i][1] = d_gamma

                    if str(gamma) == str(gamma2):
                        
                        chi = abs(a-c) + abs(b-d)
                        chi_sqr = chi**2
                        q_sqr_chi = mpq((chi_sqr - 1), chi_sqr)
                        
                        if mode == "YKR":
                            if matrix[i][i][2] == None:
                                matrix[i][i][2] = mpq(1,chi)
                            if matrix[j][j][2] == None:
                                matrix[j][j][2] = mpq(-1,chi)
                            if matrix[i][j][2] == None:
                                matrix[i][j][2] = mpq(1)
                            if matrix[j][i][2] == None:
                                matrix[j][i][2] = q_sqr_chi
                            
                        if mode == "YSR":
                                if matrix[i][i][2] == None:
                                    matrix[i][i][2] = mpq(1,chi)
                                if matrix[j][j][2] == None:
                                    matrix[j][j][2] = mpq(-1,chi)
                                if matrix[i][j][2] == None:
                                    matrix[i][j][2] = q_sqr_chi
                                if matrix[j][i][2] == None:
                                    matrix[j][i][2] = mpq(1)
                        if mode == "YOR":
                                
                                sqrt = math.sqrt(q_sqr_chi.numerator/q_sqr_chi.denominator)

                                if matrix[i][i][2] == None:
                                    matrix[i][i][2] = mpq(1,chi)
                                if matrix[j][j][2] == None:
                                    matrix[j][j][2] = mpq(-1,chi)
                                if matrix[i][j][2] == None:
                                    matrix[i][j][2] = sqrt
                                if matrix[j][i][2] == None:
                                    matrix[j][i][2] = sqrt
                    else:
                        matrix[i][j][2] = mpq(0)
                        matrix[j][i][2] = mpq(0)
                                            

        return matrix
    
    def _build_irrep(self, partition, mode="YKR"):
    

     matrix = self._generate_representation_matrix(partition, ordered_2ndlevel_tableaux(partition), mode)

     return matrix



    def _direct_sum(self, matrices):
        n = len(matrices)
        
        total_rows = sum(len(mat) for mat in matrices)
        total_cols = sum(len(mat[0]) for mat in matrices)
        
        result = [[[0, 0, 0] for _ in range(total_cols)] for _ in range(total_rows)]
        
        row_start = 0
        col_start = 0
        
        for mat in matrices:
            mat_rows = len(mat)
            mat_cols = len(mat[0])
            
            for i in range(mat_rows):
                for j in range(mat_cols):
                    result[row_start + i][col_start + j] = mat[i][j]
            
            row_start += mat_rows
            col_start += mat_cols
        
        for i in range(total_rows):
            for j in range(total_cols):
                if result[i][j] == [0, 0, 0]:
                    result[i][j] = [result[i][i][0], result[j][j][0], 0]
        
        return result

    def _build_irrep_of_transposition(self, partition, t_n: int, mode="YKR"):
        """
        Construye la representación irreducible de una transposición (como matriz por bloques).

        Args:
            partition (Snob2.IntegerPartition): La partición entera asociada.
            t_n (int): El índice de la transposición.
            mode (str, optional): El modo de la representación. Debe ser uno de 'YKR', 'YSR' o 'YOR'. Por defecto es 'YKR'.
        """
        # Caso base, la partición es de n elementos y queremos evaluar la representación en t_n
        n = sum(partition)
        if t_n == n:
            return self._build_irrep(partition, mode)
        else:
            # En otro caso, tenemos que expandir el árbol, es decir
            # Tenemos que calcular todas las transposiciones cuyos tableros de young salen de la actual
            partitions = ordered_1stlevel_tableaux(partition)
            submatrix_list = []
            for part in partitions:
                submatrix = self._build_irrep_of_transposition(part, t_n, mode)
                submatrix_list.append(submatrix)
            resp = self._direct_sum(submatrix_list)
            return resp
    
    def _calc_matrix(self, transposition):
        if self.mode == "YOR":
            self.matrices[transposition-2] = decompress_YOR(self._build_irrep_of_transposition(self.partition, transposition, self.mode))
        else:
            self.matrices[transposition-2] = decompress(self._build_irrep_of_transposition(self.partition, transposition, self.mode))

    def evaluate(self, pi):
        """
        Construye la representación irreducible de una permutación.

        Args:.
            pi (list): La permutación expresada como una lista de transposiciones adyacentes.

        Returns:
            np.ndarray: La matriz de la representación irreducible de la permutación.
        """
        transpositions = express_into_adyacent_transpositions(pi)
        
        if len(transpositions) == 0:
            return np.eye(self.matrices[0].shape[0], dtype = object) if self.mode != "YOR" else np.eye(self.matrices[0].shape[0])
        
        currMatrix = self.matrices[transpositions[0]-2]
        for transposition in transpositions[1:]:
            currMatrix = currMatrix @ self.matrices[transposition-2]

        return currMatrix
    
    def evaluate_floatingpoint(self, pi):
        """
        Construye la representación irreducible de una permutación.

        Args:.
            pi (list): La permutación expresada como una lista de transposiciones adyacentes.

        Returns:
            np.ndarray: La matriz de la representación irreducible de la permutación.
        """
        transpositions = express_into_adyacent_transpositions(pi)
        
        if len(transpositions) == 0:
            return np.eye(self.matrices[0].shape[0])
        
        currMatrix = np.array(self.matrices[transpositions[0]-2], dtype=np.float64)
        for transposition in transpositions[1:]:
            currMatrix = currMatrix @ np.array(self.matrices[transposition-2], dtype=np.float64)

        return currMatrix

    # def _mcm(self,a , b):
    #     return abs(a*b) // math.gcd(a,b)

    # def _calculate_transpositions_and_mcm(self):
    #     self.matrices = [None for _ in range(self.n - 1)]
    #     for i in range(2, self.n + 1):
    #         self._calc_matrix(i)
    #     if self.mode == "YOR":
    #         self.mcm = None
    #     else:
    #         denominators = [mat[1] for mat in self.matrices]
    #         coeficientes = [elem for matriz in denominators for fila in matriz for elem in fila]
    #         #Calculo el mcm de todos los denominadores
    #         self.mcm = reduce(self._mcm, coeficientes)
    #         #Ahora la lista denominators tiene los numeros por los que hay que multiplicar los numeradores
    #         denominators = [np.array(np.array(self.mcm / mat, dtype=int), dtype=object) for mat in denominators]
    #         numerators = [self.matrices[i][0] * denominators[i] for i in range(len(self.matrices))]
    #         #Conservamos los numeradores y el mcm
    #         self.matrices = numerators


rep = Irrep([3,2], "YKR")
for mat in rep.matrices:
    print_matrix(mat)
    print("------------------")
# pi1 = Snob2.SnElement([2,4,1,5,3,6])
# pi2 = Snob2.SnElement([3,4,5,1,2,6])
# matrix1 = rep.evaluate(pi1)
# matrix2 = rep.evaluate(pi2)
# matrix = matrix1 @ matrix2
# m = rep.evaluate(pi1*pi2)
# print(np.allclose(matrix, m, atol=1e-6))

# rho = Irrep(Snob2.IntegerPartition([3,1]), mode="YKR")
# print_matrix(rho.evaluate(Snob2.SnElement([4,1,3,2])))