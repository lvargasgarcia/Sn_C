import numpy as np
from fractions import Fraction

def compose(p1,p2):
    return [p1[p2[i]-1] for i in range(0,len(p1))]

def inverse(p):
    return [p.index(i)+1 for i in range(1,len(p) + 1)]

def williamsCondition(p,n):
    i = p.index(n)
    r = p[(i % (n-1)) + 1]
    return r % (n-1) + 1 == p[0]

def contar_decimales(num):
    # Convertimos el número a cadena, eliminamos el signo "-" si existe
    str_num = str(abs(num))
    
    # Si tiene punto decimal, contamos los dígitos después del punto
    if '.' in str_num:
        return len(str_num.split('.')[1])
    else:
        return 0

def representar_decimal(num):
    decims = contar_decimales(num)
    if decims > 0:
        numerator = int(num * 10**decims)
        denominator = int(10**decims)
        return Fraction(numerator, denominator)
    else:
        return Fraction(num)

def generate_williams_list(f, n):
    
    f_nums = []
    f_dens = []
    sequence = ""

    tau = [2,1] + [i for i in range(3,n+1)] # transposición (1,2)
    sigma = [(i+1) for i in range(1,n)] + [1] # Ciclo completo
    q = [(n-i+1) for i in range(1,n+1)] # El 1 va al n, el 2 al n-1, etc.
    qtau = compose(q, tau)
    qsigma = compose(q, sigma)
    qsigmatau = compose(qsigma, tau)
    invSigma = inverse(sigma)
    
    p = compose(qsigma, tau)

    f_nums.append(representar_decimal(f(tuple(p))).numerator)
    f_dens.append(representar_decimal(f(tuple(p))).denominator)


    while(p != qtau):
        
        if(p != qsigmatau):
            if(williamsCondition(p,n) and p != qsigma):
                p = compose(p, tau)
                sequence += "t"
            else:
                p = compose(p, invSigma)
                sequence += "i"
        else:
            p = compose(p, invSigma)
            sequence += "i"
        
        f_nums.append(representar_decimal(f(tuple(p))).numerator)
        f_dens.append(representar_decimal(f(tuple(p))).denominator)

    return (np.array(f_nums), np.array(f_dens), sequence)
