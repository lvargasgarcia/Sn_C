from sympy.combinatorics import Permutation


# En este fichero se implementa la descomposición de una permutación de S_n en producto de transposiciones adyacentes
# La idea es tener una funcion con esta signatura express_into_adyacent_transpositions :: SnElement -> [int], donde int
# Es la lista ordenada de enteros que representan los indices de las transposiciones
#                      (1 2 3 4 5)
# Por ejemplo, si pi = ----------- entonces express_into_adyacent_transpositions(pi) = [3,5] porque pi = (2 3) (4 5)
#                      (1 3 2 5 4)

# def remove_consecutive_duplicates(lst):
#     if not lst:
#         return []
    
#     result = [lst[0]]
#     for i in range(1, len(lst)):
#         if lst[i] != lst[i - 1]:
#             result.append(lst[i])
    
#     return result

def get_list_from_SnElement(pi):
    return list(map(int, str(pi).strip("[]").split()))

# def image(pi, n):
#     return pi[n-1]

# def express_into_disjoint_cycles(pi):
#     pi = get_list_from_SnElement(pi)
#     n = len(pi)
#     visited_dict = {i: False for i in range(1, n+1)}
#     cycles = []
#     for i in range(1, n+1):
#         # Si el número no ha sido visitado, entonces comenzamos un nuevo ciclo
#         if not visited_dict[i]:
#             cycle = []
#             j = i
#             while not visited_dict[j]:
#                 visited_dict[j] = True
#                 cycle.append(j)
#                 j = image(pi, j)
#             cycles.append(cycle)
#     return cycles

# def express_cycle_into_transpositions(cycle):
#     if len(cycle) == 1:
#         return None
#     transpositions = []
#     initialElem = cycle[0]
#     for elem in reversed(cycle):
#         if elem == cycle[0]:
#             break
#         transpositions.append((initialElem, elem))
#     return transpositions

# def express_into_transpositions(pi):
#     cycles = express_into_disjoint_cycles(pi)
#     transpositions = []
#     for cycle in cycles:
#         trans = express_cycle_into_transpositions(cycle)
#         if trans != None:
#             transpositions += trans
#     return transpositions

# Esta función sirve para limpiar una lista de transposiciones adyacentes, ya que una transposicion es su propia inversa

def reduce(list):
    i = 0
    reduced = False
    while i < len(list)-1:
        if list[i] == list[i+1]:
            list.pop(i)
            list.pop(i)
            reduced = True
        i = i + 1
    return reduced
    

def clean(list):
    return clean(list) if reduce(list) else list


def express_transposition_into_adyacent_transpositions_right(transposition):
    (a,b) = transposition
    if a + 1 == b:
        return [b]
    else:
        return express_transposition_into_adyacent_transpositions_right((a+1,b)) + [a+1]
    
def express_transposition_into_adyacent_transpositions(transposition):
    list_right = express_transposition_into_adyacent_transpositions_right(transposition)
    list_left = list_right[1:][::-1]
    return list_left + list_right

def express_into_adyacent_transpositions(pi):
    transpositions = Permutation([0] + pi).transpositions()
    adyacent_transpositions = []
    for transposition in transpositions:
        adyacent_transpositions += express_transposition_into_adyacent_transpositions(transposition)
    print(adyacent_transpositions)
    return clean(adyacent_transpositions)
    
print(express_into_adyacent_transpositions([10,2,1,4,3,5,6,7,8,9]))

# # Algunas pruebas:
# n = 4
# G = Snob2.Sn(n)
# for i in range(len(G)):
#     print("--------------------")
#     pi = [0] + get_list_from_SnElement(G[i])
#     print(pi[1:])
#     transpositions = Permutation(pi).transpositions()
#     print(f"La permutación {pi} se descompone en las transposiciones {transpositions}")
#     adyacent_transpositions = express_into_adyacent_transpositions(pi)
#     print(f"La permutación {pi} se descompone en las transposiciones adyacentes {adyacent_transpositions}")
#     print("--------------------")

# pi = Snob2.SnElement([2,3,4,1])
# adyacent_transpositions = express_into_adyacent_transpositions(pi)
# print(f"La permutación {pi} se descompone en las transposiciones adyacentes {adyacent_transpositions}")
# print("--------------------")




