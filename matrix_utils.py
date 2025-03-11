import numpy as np


def apply_correction_transformation(matriz):
        n = len(matriz)
        return [[matriz[n - 1 - j][n - 1 - i] for j in range(n)] for i in range(n)]

def decompress_YOR(matrix):
    """
    Convierte una matriz de fracciones por bloques en una matriz de números reales.

    Args:
        matrix (list): La matriz de fracciones por bloques a convertir.

    Returns:
        np.ndarray: La matriz convertida a números reales.
    """
    # Inicializamos una lista vacía para los bloques
    blocks = []
    
    for row in matrix:
        block_row = []
        for block in row:
            # Obtenemos las dimensiones y el coeficiente
            block_rows, block_cols, coef = block

            if block_rows == block_cols:
                # Si el bloque es cuadrado, inicializamos como la matriz identidad
                block_matrix = np.eye(block_rows) * coef
            else:
                block_matrix = np.zeros((block_rows, block_cols))
            
            block_row.append(block_matrix)
        
        blocks.append(np.hstack(block_row))
    result_matrix = np.vstack(blocks)
    
    return np.array(apply_correction_transformation(result_matrix))

def decompress(matrix):
    """
    Convierte una matriz de fracciones por bloques en una matriz de números reales.

    Args:
        matrix (list): La matriz de fracciones por bloques a convertir.

    Returns:
        np.ndarray: La matriz convertida a números reales.
    """
    # Inicializamos una lista vacía para los bloques
    blocks = []
    
    for row in matrix:
        block_row = []
        for block in row:
            # Obtenemos las dimensiones y el coeficiente
            block_rows, block_cols, coef = block

            if block_rows == block_cols:
                # Si el bloque es cuadrado, inicializamos como la matriz identidad
                block_matrix = np.eye(block_rows, dtype=object) * coef
            else:
                block_matrix = np.zeros((block_rows, block_cols), dtype=object)
            
            block_row.append(block_matrix)
        
        blocks.append(np.hstack(block_row))
    result_matrix = np.vstack(blocks)
    
    return result_matrix

# def decompress(matrix):
#     """
#     Convierte una matriz de fracciones por bloques en una matriz de números reales.
#     La diferencia con el modo YOR es que aquí devolvemos tuplas con matrices(para numeradores y denominadores).

#     Args:
#         matrix (list): La matriz de fracciones por bloques a convertir.

#     Returns:
#         np.ndarray: La matriz convertida a números reales.
#     """
#     # Inicializamos una lista vacía para los bloques
#     nums_blocks = []
#     dens_blocks = []
    
#     for row in matrix:
#         num_block_row = []
#         den_block_row = []
#         for block in row:
#             # Obtenemos las dimensiones y el coeficiente
#             block_rows, block_cols, coef = block

#             if block_rows == block_cols:
#                 # Si el bloque es cuadrado, inicializamos como la matriz identidad
#                 num_block_matrix = np.eye(block_rows) * coef.numerator
#                 den_block_matrix = np.ones((block_rows, block_cols)) * coef.denominator
#             else:
#                 num_block_matrix = np.zeros((block_rows, block_cols))
#                 den_block_matrix = np.ones((block_rows, block_cols))
            
#             num_block_row.append(num_block_matrix)
#             den_block_row.append(den_block_matrix)
        
#         nums_blocks.append(np.hstack(num_block_row))
#         dens_blocks.append(np.hstack(den_block_row))
#     result_tuple = (np.vstack(nums_blocks).astype(int), np.vstack(dens_blocks).astype(int))
    
#     return result_tuple

def pretty_print(matrix):
    """
    Imprime una matriz de representación por bloques de forma legible.

    Args:
        matrix (list): La matriz de representación por bloques a imprimir.
    """
    print("Matriz de representación:")
    
    print("---------Dimensiones---------")
    
    dimensions = [[str((cell[0], cell[1])) for cell in row] for row in matrix]
    col_widths = [max(len(row[i]) for row in dimensions) for i in range(len(dimensions[0]))]
    
    for row in dimensions:
        print("  ".join(cell.rjust(col_widths[i]) for i, cell in enumerate(row)))
    
    print("---------Coeficientes---------")

    coefficients = [[str(cell[2]) for cell in row] for row in matrix]
    col_widths = [max(len(row[i]) for row in coefficients) for i in range(len(coefficients[0]))]
    
    for row in coefficients:
        print("  ".join(cell.rjust(col_widths[i]) for i, cell in enumerate(row)))

def print_matrix(matrix):
    """
    Imprime una matriz de números reales de forma legible.

    Args:
        matrix (list): La matriz de números reales a imprimir.
    """
    for sublista in matrix:
        print(f"[{' '.join(map(str, sublista))}]")

