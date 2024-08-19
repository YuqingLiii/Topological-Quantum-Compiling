import sympy as sp
import numpy as np
     
def UnitaryMatrix(list):
    assert len(list)==4 or len(list)==16 , "please check the input matrix"
            
    if len(list)==4:
        matrix=sp.Matrix(np.array(list).reshape(2, 2))
    else :
        matrix=sp.Matrix(np.array(list).reshape(4, 4))
    
    assert abs(matrix.det())==1, "input matirx is not unitary"
    return matrix