import sympy as sp
from UnitaryMatrix import UnitaryMatrix

def gate_distance(guess_coeff, target_matrix):
    a, b, c, d = guess_coeff
    tau = (sp.sqrt(5) - 1) / 2
    F = sp.Matrix([
        [tau, sp.sqrt(tau)],
        [sp.sqrt(tau), -tau]
    ])
    
    estim_matrix = sp.exp(sp.I * a) * Rz(b) * F * Rz(c) * F * Rz(d)
    
    difference_matrix = estim_matrix - target_matrix
    frobenius_norm = difference_matrix[0,0]**2+difference_matrix[0,1]**2+difference_matrix[1,0]**2+difference_matrix[1,1]**2
    
    return frobenius_norm.evalf()

def gate_distance_single(estim_matrix, target_matrix):
    
    difference_matrix = estim_matrix - target_matrix
    frobenius_norm = difference_matrix[0,0]**2+difference_matrix[0,1]**2+difference_matrix[1,0]**2+difference_matrix[1,1]**2
    
    return frobenius_norm.evalf()

def Rz(angle):
        return sp.Matrix([[sp.exp(-0.5*angle*sp.I),0],[0,sp.exp(0.5*angle*sp.I)]])

def matrix_decomposition(target_matrix):
    tau=(sp.sqrt(5)-1)/2
    omega=sp.exp(sp.pi*sp.I/5)
    F_matrix=UnitaryMatrix([tau,sp.sqrt(tau),sp.sqrt(tau),-tau])
    #decomposition_form=sp.exp(delta*sp.I)*Rz(alpha)*F_matrix*Rz(beta)*F_matrix*Rz(gamma)
    if sp.im(target_matrix[0,0])==0 and sp.im(target_matrix[1,1])==0:
        if sp.re(target_matrix[0,0])==sp.re(target_matrix[1,1]):
            delta=0
        elif sp.re(target_matrix[0,0])==-sp.re(target_matrix[1,1]):
            delta=sp.pi/2
        else:
            print("不能分解这个门")
    elif sp.im(target_matrix[0,0])-sp.im(target_matrix[1,1])==0:
        tan_delta=-(sp.im(target_matrix[0,0])+sp.im(target_matrix[1,1]))/(sp.re(target_matrix[0,0])+sp.re(target_matrix[1,1]))
        delta=sp.atan(tan_delta)
    else:
        tan_delta=(sp.re(target_matrix[0,0])-sp.re(target_matrix[1,1]))/(sp.im(target_matrix[0,0])-sp.im(target_matrix[1,1]))
        delta=sp.atan(tan_delta)
  
    target_matrix_=sp.exp(-sp.I*delta)*target_matrix

    assert target_matrix_[0,0].conjugate()==target_matrix_[1,1],"相位提取失败"
    
    flag=0
    if(target_matrix[0,0]!=0):
        
        cos_beta=(target_matrix_[0,0]*target_matrix_[1,1]-tau**4-tau**2)/(2*tau**3)
        beta=sp.acos(cos_beta).evalf()
        
        exp_a_sub_c=target_matrix_[1,0]*tau**(-3/2)/(-2*sp.I)/sp.sin(beta/2)      
        
        a_sub_c = 2*sp.atan2(sp.im(exp_a_sub_c), sp.re(exp_a_sub_c))
        exp_a_add_c=target_matrix_[0,1]/(sp.exp(sp.I*beta/2)*tau+sp.exp(-sp.I*beta/2)*tau**2)
        a_add_c = -2*sp.atan2(sp.im(exp_a_add_c), sp.re(exp_a_add_c))
        alpha=(a_sub_c+a_add_c)/2
        gamma=a_add_c-alpha
        result=[delta.evalf(),alpha.evalf(),beta.evalf(),gamma.evalf()]
        return result,flag
    else:
        flag=1
        
    return result,flag

print(1)
H=UnitaryMatrix([1/sp.sqrt(2),1/sp.sqrt(2),1/sp.sqrt(2),-1/sp.sqrt(2)])
print(gate_distance([sp.pi/4,1.81618578660117,1.62984763439859,1.81618578660117],H))