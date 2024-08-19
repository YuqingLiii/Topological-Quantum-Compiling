from TwoRings import RingCycIntegers,SubRealRing
from UnitaryMatrix import UnitaryMatrix
import cmath
import sympy as sp
import numpy as np
import random
from GateDecomposion import matrix_decomposition,gate_distance
from wolframclient.language import wl, wlexpr
from wolframclient.evaluation import WolframLanguageSession
from timeit import default_timer as timer

def gate_distance_single(estim_matrix, target_matrix):
    
    difference_matrix = estim_matrix - target_matrix
    frobenius_norm = difference_matrix[0,0]*difference_matrix[0,0].conjugate()+difference_matrix[0,1]*difference_matrix[0,1].conjugate()+difference_matrix[1,0]*difference_matrix[1,0].conjugate()+difference_matrix[1,1]*difference_matrix[1,1].conjugate()
    
    return frobenius_norm.evalf()
class NTA:
    
    tau=(sp.sqrt(5)-1)/2
    omega=sp.exp(sp.pi*sp.I/5)
    
    F_matrix=UnitaryMatrix([tau,sp.sqrt(tau),sp.sqrt(tau),-tau])
    T_matrix=UnitaryMatrix([1,0,0,omega])
    
    pauli_1=UnitaryMatrix([omega**6,0,0,omega**13])
    pauli_2=F_matrix*pauli_1*F_matrix
    
    def __init__(self,target_matrix,precision):
        self.target_matrix=UnitaryMatrix(target_matrix)
        self.precision=precision
       
    def Rz(self,angle):
        return sp.Matrix([[sp.exp(-0.5*angle*sp.I),0],[0,sp.exp(0.5*angle*sp.I)]])
    
    def fibonacci(self,n):
        a, b = 0, 1
        for _ in range(n):
            a, b = b, a + b
        return a

    def approx_real(self,x, n):
        
        p = self.fibonacci(n)
        q = self.fibonacci(n + 1)
        u = (-1)**(n + 1) * p
        v = (-1)**n * self.fibonacci(n - 1)
        c = int(x * q)
        a = c * v + p * int(c * u / q)
        b = c * u - q * int(c * u / q)
        #print(a,b)
        return SubRealRing([a,b])
    
    def random_sample(self,theta,epsilon,r):
        
        phi=1/self.__class__.tau
        C=sp.sqrt(phi/(4*r))
        m=sp.ceiling(sp.log(C*epsilon*r,self.__class__.tau))+1
        N=sp.ceiling(phi**m)
        
        y_min = r * phi**m * (sp.sin(theta) - epsilon * (sp.sqrt(4 - epsilon**2) * sp.cos(theta) + epsilon * sp.sin(theta)) / 2)
        y_max = r * phi**m * (sp.sin(theta) + epsilon * (sp.sqrt(4 - epsilon**2) * sp.cos(theta) - epsilon * sp.sin(theta)) / 2)
        x_max = r * phi**m * ((1 - epsilon**2 / 2) * sp.cos(theta) - epsilon * sp.sqrt(1 - epsilon**2 / 4) * sp.sin(theta))
        x_c = x_max - r * epsilon**2 * phi**m / (4 * sp.cos(theta))

        j = random.randint(1, N-1)
        
        y = y_min + j * (y_max - y_min) / N
 
        a_y_tau_b_y = self.approx_real(y / sp.sqrt(2 - self.__class__.tau), m)
        a_y_tau_b_y_ring=a_y_tau_b_y.transform_ringcycinteger()
        
        x = x_c - ((a_y_tau_b_y).get_value() * sp.sqrt(2 - self.__class__.tau) - y_min) * sp.tan(theta)
        a_x_tau_b_x = self.approx_real(x, m)
        a_x_tau_b_x_ring=a_x_tau_b_x.transform_ringcycinteger()

        temp=RingCycIntegers([-1,2,-1,1])
        u0=temp.mul(a_y_tau_b_y_ring).add(a_x_tau_b_x_ring)
        
        return u0
            
    def extract_a_b(self,xi):
        
        a, b = sp.symbols('a, b')
        
        equation = sp.Eq(xi, a + b * self.__class__.tau)
        solutions = sp.solve(equation, (a, b))
        
        return solutions.get(a), solutions.get(b)

    def easy_solvable(self,fl):
        for i in range(0,len(fl)):
            xi=fl[i][0]
            k=fl[i][1]
            
            if k%2 == 1:
                if xi.get_value()!=5:
                    
                    p=xi.N_map()            
                    r=int(p.get_value()%5)
                    #print(p.get_value())
                    #print(r!=1)
                    if (sp.isprime(int(p.get_value()))==False) or (r!=0 and r!=1):
                        
                        return False
        return True
             
    def easy_factor(self,xi):
        
        a, b = xi.coeffs[0],xi.coeffs[1]
        
        c = sp.gcd(a, b)
        a1 = a // c
        b1 = b // c

        xi1 = SubRealRing([a1,b1])
        
        d = int(sp.sqrt(c))
        if d * d == c:
            ret = [(SubRealRing([d,0]), 2)]
        else:
            d = int(sp.sqrt(c / 5))
            if d * d * 5 == c:
                ret = [(SubRealRing([d,0]), 2), (SubRealRing([5,0]), 1)]
            else:
                return [(xi, 1)]
        
        n=xi1.N_map()
        
        if n.get_value() % 5 == 0:
            xi2_a=(3*xi1.get_coeff(0)+xi1.get_coeff(1))//5
            xi2_b=(xi1.get_coeff(0)+2*xi1.get_coeff(1))//5
            
            assert xi2_a%1==0 and xi2_b%1==0,"系数不是整数"
            xi2 = SubRealRing([xi2_a,xi2_b])

            ret.append((SubRealRing([2,-1]), 1))
            ret.append((xi2,1))

            return ret
        else:
            ret.append((xi1, 1))
            
            return ret
    
    def automorphism(xi):
        a,b=extract_a_b(xi)
        return a-b*(self.__class__.tau+1)
    
    def tonelli_shanks(self,n,p):
        #assert pow(n, (p - 1) // 2, p) == 1, "n不是模数p的二次剩余"
    
        if n == 0:
            return 0
        if p == 2:
            return n

        s = 0
        q = p - 1
        while q % 2 == 0:
            s += 1
            q //= 2
        
        if s == 1:
            
            r = pow(n, int((p + 1) // 4), int(p))
            return r, p - r
        
        z = 2
        
        while pow(z, int((p - 1) // 2), int(p)) !=  int(p)-1:
             z += 1
       
        c = pow(z, int(q), int(p))
        r = pow(n, int((q + 1) // 2), int(p))
        t = pow(n, int(q), int(p))
        m = s
        
        while t != 1:
            
            t2i = t
            i = 0
            for i in range(1, m):
                t2i = pow(t2i, 2, int(p))
                if t2i == 1:
                    break
            b = pow(c, 2 ** (m - i - 1), int(p))
            r = (r * b) 
            t = (t * b * b) 
            c = (b * b) 
            m = i
        
        return r,p-r
    
    def splitting_root(self,xi):
        a,b=xi.get_coeff(0),xi.get_coeff(1)
        p=xi.N_map()
        assert b%p.get_value()!=0 , "b不能是模数p的倍数"
        
        b1=sp.mod_inverse(b, int(p.get_value()))

        return self.tonelli_shanks(-a*b1-2,p.get_value())
    
    def is_Divide(self,item):
        ua=item.get_coeff(0)-item.get_coeff(1)+item.get_coeff(2)-item.get_coeff(3)
        return ua%5==0
        #re_matrix=sp.Matrix([[4, 1, -1, 1], [-3, 3, 2, -2], [2, -2, 2, 3], [-1, 1, -1, 1]])
        #solution=re_matrix*sp.Matrix([[item.get_coeff(0)],[item.get_coeff(1)],[item.get_coeff(2)],[item.get_coeff(3)]])
        
        #is_integer_matrix = all(element%5==0 for element in solution)
        
        #return is_integer_matrix
        
    def factor_1_plus_omega(self,item):
        re_matrix=sp.Matrix([[4, 1, -1, 1], [-3, 3, 2, -2], [2, -2, 2, 3], [-1, 1, -1, 1]])
        solution=re_matrix*sp.Matrix([[item.get_coeff(0)],[item.get_coeff(1)],[item.get_coeff(2)],[item.get_coeff(3)]])
        
        result=RingCycIntegers([int(solution[0,0]//5),int(solution[1,0]//5),int(solution[2,0]//5),int(solution[3,0]//5)])
       
        return result
    
    def mod_inverse(self,xi):

        a=xi.get_coeff(0)%5
        b=xi.get_coeff(1)%5
        c=xi.get_coeff(2)%5
        d=xi.get_coeff(3)%5

        with WolframLanguageSession("D:\\mathematica\\WolframKernel.exe") as session:
            result = session.evaluate(wlexpr( "FindInstance[{}*m + {}*n + {}*x + {}*y == 1, {{x, y, m, n}}, Integers]".format((a+c-b-d)%5, (-a-c+b-4*d)%5, (a-4*c-b-d)%5, (-a-c-4*b+4*d)%5)))
            print((a+c-b-d)%5, (-a-c+b-4*d)%5, (a-4*c-b-d)%5, (-a-c-4*b+4*d)%5)
            """
            x=self.custom_mod5(result[0][0][1])
            y=self.custom_mod5(result[0][1][1])
            m=self.custom_mod5(result[0][2][1])
            n=self.custom_mod5(result[0][3][1])
            """
            x=result[0][0][1]%5
            y=result[0][1][1]%5
            m=result[0][2][1]%5
            n=result[0][3][1]%5
        print([m,n,x,y])
        #print(self.is_Divide(RingCycIntegers([m,n,x,y]).mul(xi).sub(RingCycIntegers([1,0,0,0]))))
        return RingCycIntegers([m,n,x,y])
             
    def custom_mod5(self,n):
        remainder = n % 5
        if remainder == 3:
            return -2
        elif remainder == 4:
            return -1
        else:
            return remainder
              
    def binary_gcd(self, a, b, log_file="gcd_log.txt"):
    # 打开文件进行写操作，如果文件不存在则会创建
        
        with open(log_file, "a", encoding="utf-8") as f:
            def log_message(message):
                f.write(message +"\n")
            
            log_message("本次的a，b是：")
            log_message(a.get_str())
            log_message(b.get_str())
        
            """
                在 Z[ω] 环上实现 BINARY-GCD 算法
            """
            #print("本次的a，b是：")
            #print(a.get_str())
            #print(b.get_str())
            if a.get_value() == 0:
                return b
            if b.get_value() == 0:
                return a
            log_message(f"is_Divide(a): {self.is_Divide(a)}, is_Divide(b): {self.is_Divide(b)}")
            #print(f"is_Divide(a): {self.is_Divide(a)}, is_Divide(b): {self.is_Divide(b)}")
            if self.is_Divide(a) == True and self.is_Divide(b) == True:
                a1 = self.factor_1_plus_omega(a)
                b1 = self.factor_1_plus_omega(b)
                return RingCycIntegers([1,1,0,0]).mul(self.binary_gcd(a1, b1, log_file))

            if self.is_Divide(a) == True and self.is_Divide(b) == False:
                a1 = self.factor_1_plus_omega(a)
                log_message(f"a可以分且a1是{a1.get_str()}\n")
                return self.binary_gcd(a1,b, log_file)
            
            if self.is_Divide(a) == False and self.is_Divide(b) == True:
                b = self.factor_1_plus_omega(b)
                
                while self.is_Divide(b) == True:
                    b = self.factor_1_plus_omega(b)
                log_message(f"b可以分且b是{b.get_str()} \n")
                return self.binary_gcd(a,b, log_file)
                
            u = v = RingCycIntegers([1,0,0,0])
            if self.is_Divide(a) == False and self.is_Divide(b) == False:
                
                a_mod = a.get_coeff(0) - a.get_coeff(1) + a.get_coeff(2) - a.get_coeff(3)
                if a_mod % 5 == 2:
                    u = RingCycIntegers([-1, 1, 0, 0])
                elif a_mod % 5 == 3:
                    u = RingCycIntegers([1, -1, 0, 0])
                else:
                    u = self.custom_mod5(sp.mod_inverse(a_mod, 5))
                    u = RingCycIntegers([u, 0, 0, 0])
                    
                #u = self.custom_mod5(sp.mod_inverse(a_mod, 5))
                #u = RingCycIntegers([u, 0, 0, 0])
                b_mod = b.get_coeff(0) - b.get_coeff(1) + b.get_coeff(2) - b.get_coeff(3)
                if b_mod % 5 == 2:
                    v = RingCycIntegers([-1, 1, 0, 0])
                elif b_mod % 5 == 3:
                    v = RingCycIntegers([1, -1, 0, 0])
                else:
                    v =self.custom_mod5(sp.mod_inverse(b_mod, 5))
                    v = RingCycIntegers([v, 0, 0, 0])
                #v =self.custom_mod5(sp.mod_inverse(b_mod, 5))
                #v = RingCycIntegers([v, 0, 0, 0])
                #u=self.mod_inverse(a)
                #v=self.mod_inverse(b)
                log_message("计算出来的逆元都是：")
                log_message(f"{u.get_str()}, {v.get_str()}")
                #print("计算出来的逆元都是：")
                #print(f"{u.get_str()}, {v.get_str()}")
            log_message("a和b的模大小是：")
            log_message("   ")
            log_message(str(a.abs_quadratic().get_value()))
            log_message("   ")
            log_message(str(b.abs_quadratic().get_value()))
            log_message(" \n  ")
            if a.abs_quadratic().get_value() < b.abs_quadratic().get_value():
                c = a
            else:
                c = b

            return self.binary_gcd(c, u.mul(a).sub(v.mul(b)), log_file)

    def unit_dlog(self,u):
        """实现 UNIT-DLOG 算法"""
        
        a, b = u.get_coeff(0),u.get_coeff(1)
        
        s = 1
        k = 0
        if a < 0:
            a = -a
            b = -b
            s = -s
        mu = a * b
        while abs(mu) > 1:
            #print(a,b,mu)
            with open("unit_dlog", "a", encoding="utf-8") as f:
                f.write(str(a))
                f.write("   ")
                f.write(str(b))
                f.write("   ")
                f.write(str(mu))
                f.write(" \n  ")
            temp=a-b
            if mu>1:
                a=b
                b=temp
                k -= 1
            else:
                #a, b = a, temp
                a=a
                b=temp
                k += 1
            mu = a * b
        print("place 1")
        print(mu)
        tau=(sp.sqrt(5)-1)/2
        v_values = [1, -1, tau, -tau, tau**2, -tau**2]
        for v in v_values:
            if sp.simplify(u.get_value() - s * v * tau**k) < 10**(-2):                    
                if v == -1:
                    s = -s
                elif v == tau:
                    k += 1
                elif v == -tau:
                    s = -s
                    k += 1
                elif v == tau**2:
                    k += 2
                elif v == -tau**2:
                    s = -s
                    k += 2
                return s, k

    def solve_norm_equation(self,xi):
        xi_dot= xi.automorphism()
        if xi.get_value()<0 or xi_dot.get_value()<0:
            return False
        else:
            factor_list=self.easy_factor(xi)
        if not self.easy_solvable(factor_list):
            return False
        x=RingCycIntegers([1,0,0,0])
        for i in range(0,len(factor_list)):
            print(i)
            xi_i=factor_list[i][0]
            #xi_i=SubRealRing([15,-8])
            m=factor_list[i][1]
            x=x.mul(xi_i.pow(sp.floor(m/2)).transform_ringcycinteger())
            print(xi_i.get_str())
            if m%2==1:
                if xi_i.get_value()==5:
                    x=x.mul(SubRealRing([1,2]).transform_ringcycinteger())
                else:
                    if xi_i.get_value()==SubRealRing([2,-1]).get_value():
                        x=x.mul(RingCycIntegers([-1,2,-1,1]))
                    else:
                        print("%%")
                        M,_=self.splitting_root(xi_i)
                        print(M)
                        #y=self.binary_gcd(xi_i.transform_ringcycinteger(),RingCycIntegers([M+1,-2,1,-1]))
                        
                        y=RingCycIntegers([1,0,0,0])
                        y_qua=y.abs_quadratic(True)
                        a=xi_i.get_coeff(0)
                        b=xi_i.get_coeff(1)
                        c=y_qua.get_coeff(0)
                        d=y_qua.get_coeff(1)
                        n=int((b*c-a*d)/(c**2-d**2-c*d))
                        m=int((a-n*d)/c)
                        print(m,n)
                        u=SubRealRing([m,n])
                        s,m=self.unit_dlog(u)
                        x=x.mul(SubRealRing([0,1]).transform_ringcycinteger.pow(m//2)).mul(y)
        return x
           
    def gauss_complexity(self,u):
        """我的脑子有问题了，这个函数不用看了"""
        if u==0:
            return 0
        max_k=10
        for k in range(max_k):
            if sp.simplify(self.__class__.omega**k - u) == 0:
                return 2
            if sp.simplify(self.__class__.tau*self.__class__.omega^k-u)==0 or sp.simplify(self.__class__.tau**(-1)*self.__class__.omega^k-u)==0  :
                return 3
        omega=self.__class__.omega
        expanded_mu = sp.expand_complex(u)
        collected_mu = sp.collect(expanded_mu, omega)

        b = collected_mu.coeff(omega**1)
        c = collected_mu.coeff(omega**2)
        d = collected_mu.coeff(omega**3)
        a=u-b*omega + c*omega**2 + d*omega**3
        max_G=2.5*(a**2+b**2+c**2+d**2)
        assert max_G>=5, "高斯复杂度计算有问题"
        return sp.ceiling(max_G)
           
    def FTV(self,j,V):
        omega_j=RingCycIntegers([0,1,0,0]).pow(j)
        tau_ring=SubRealRing([0,1]).transform_ringcycinteger()
        a=(V[0][0].add(V[1][0].mul(omega_j))).mul(tau_ring)
        c=V[0][0].sub(V[1][0].mul(omega_j.mul(tau_ring)))
        b=V[0][1].mul(tau_ring).add(omega_j.mul(V[1][1]))
        d=(V[0][1].sub(omega_j.mul(V[1][1]))).mul(tau_ring)
        return [[a,b],[c,d]]
    
    def exact_synthesize(self,u, v, k):
        g=u.automorphism().abs_quadratic().get_value()
        tau_ring=SubRealRing([0,1]).transform_ringcycinteger()
        C=""
        C_Matrix=sp.Matrix([[1,0],[0,1]])
        #print(u.get_str())
        #print(v.get_str())
        V=[[u,v.conjugate_()],[v,u.minus().conjugate_()]]
        J_list=[]
      
        while g>=2:
            #print(g)
            J = min(range(1, 11), key=lambda j: self.FTV(j,V)[0][0].automorphism().abs_quadratic().get_value())

            V=self.FTV(J,V)
            
            g=V[0][0].automorphism().abs_quadratic().get_value()
            C=f"F T^{10-J} "+C
            C_Matrix=self.__class__.F_matrix*self.__class__.T_matrix**(10-J)*C_Matrix

        flag=V[0][0]
        if flag.get_coeff(0)==-1 and flag.get_coeff(1)==0:
            k=5
        elif flag.get_coeff(0)==1 and flag.get_coeff(1)==0:
            k=0
        elif flag.get_coeff(1)==1 and flag.get_coeff(0)==0:
            k=1
        elif flag.get_coeff(1)==-1 and flag.get_coeff(0)==0:
            k=6
        elif flag.get_coeff(2)==-1 and flag.get_coeff(1)==0:
            k=7
        elif flag.get_coeff(2)==1 and flag.get_coeff(1)==0:
            k=2
        elif flag.get_coeff(3)==1 and flag.get_coeff(1)==0:
            k=3
        elif flag.get_coeff(3)==-1 and flag.get_coeff(1)==0:
            k=8
        else:
            if flag.get_coeff(0)==-1:
                k=4
            else:
                k=9
        d_pie=sp.arg(u.minus().conjugate_().get_value())
        d=sp.arg(C_Matrix[1,1].evalf())
        j=int((d_pie-d)/(0.628318530717959))%9

        C=f"omega^{k} T^{j} "+C
        C_Matrix=self.__class__.omega**(k)*self.__class__.T_matrix**(j)*C_Matrix
        C_Matrix_=sp.Matrix([[C_Matrix[0,0].evalf(),C_Matrix[0,1].evalf()],[C_Matrix[1,0].evalf(),C_Matrix[1,1].evalf()]])
        return C,C_Matrix_
    
    def QuantumCompile_1qubit(self,phi):
        C=sp.sqrt(self.__class__.tau**(-1)/4)
        epsilon=self.precision
        m=sp.ceiling(sp.log(C*epsilon,self.__class__.tau))+1
        
        theta=0
        k=0
        #print(phi)
        for i in range(-10,10):
            theta=-phi/ 4 - (sp.pi / 5) *i
            if theta<=sp.pi/5 and theta>=0:
                k=i
                break
        
        not_found=True
        u=0
        v=0
   
        while not_found:
            #u_0=RingCycIntegers([22,-10,28,-14])
            u_0=self.random_sample(theta,epsilon,1)
            temp_phi=SubRealRing([1,1])
            xi = temp_phi.mul(temp_phi.pow(2*m).sub(u_0.abs_quadratic(True)))
            #print(sp.sqrt(u_0.abs_quadratic().get_value()))
            #print(temp_phi.pow(m).get_value())
            #print(u_0.get_str())
            factor_list = self.easy_factor(xi)
            #not_found = False
            #for item in factor_list:
            #    print(item[0].get_str())
            if self.easy_solvable(factor_list):
                print("抽样结束, u0: {}".format(u_0.get_str()))
                #print(u_0.get_str())
                #print(temp_phi.pow(m).get_value())
                not_found = False
                temp_tau=SubRealRing([0,1]).transform_ringcycinteger()
                temp_omega=RingCycIntegers([0,1,0,0])
                u = temp_omega.pow(k).mul(temp_tau.pow(m)).mul(u_0) 
                #print(temp_omega.pow(k).get_str())
                #print(temp_tau.pow(m).get_value())
                #print(temp_tau.pow(m).get_str())
                #print((temp_tau.pow(m)).abs_quadratic().get_value())
                #print(u.abs_quadratic().get_value())
                a=xi.get_coeff(0)
                b=xi.get_coeff(1)
                with WolframLanguageSession("D:\\mathematica\\WolframKernel.exe") as session:
                    
                    result = session.evaluate(wlexpr( "Solve[{{x^2 + y^2 + m^2 + n^2 + x*y + y*m + m*n == {},  x*y + x*m - x*n + y*m + y*n + m*n == {}}}, {{x, y, m, n}}, Integers]".format(a,b)))
                    x=result[0][0][1]
                    y=result[0][1][1]
                    m_=result[0][2][1]
                    n=result[0][3][1]
                solution=RingCycIntegers([x,y,m_,n])
                v = temp_tau.pow(m).mul(solution)

                print("解方程结束")
        C,C_Matrix = self.exact_synthesize(u, v, 0)

        return C,C_Matrix
            
    def quantum_complie(self):
        parameters,flag=matrix_decomposition(self.target_matrix)
        
        FT_circuit=""
        if flag==0:
            
            alpha=parameters[1]
            beta=parameters[2]
            gamma=parameters[3]
            print("="*86)
            print(f"***********开始分解第一个门，相位是 {alpha}*********************************")
            time1=timer()
            FT_circuit,decompos_gate=self.QuantumCompile_1qubit(alpha)
            
            FT_circuit=FT_circuit+"F"
            decompos_gate=decompos_gate*self.__class__.F_matrix
            time2=timer()
            print(f"第一个门分解完成，用时 {time2-time1:.4f} 秒, 门序列是 {FT_circuit}")

            # 开始分解第二个门
            print("="*86)
            print(f"***********开始分解第二个门，相位是 {beta}*********************************")
            FT_circuit_sec,gate=self.QuantumCompile_1qubit(beta)
            FT_circuit=FT_circuit+FT_circuit_sec+"F"
            decompos_gate=decompos_gate*gate*self.__class__.F_matrix
            
            time3=timer()
            print(f"第二个门分解完成，用时 {time3-time2:.4f} 秒, 门序列是 {FT_circuit}")
            
            # 开始分解第三个门
            print("="*86)
            print(f"***********开始分解第三个门，相位是 {gamma}*********************************")
            FT_circuit_third,gate2=self.QuantumCompile_1qubit(gamma)
            FT_circuit=FT_circuit+FT_circuit_third
            decompos_gate=decompos_gate*gate2
            
            time4=timer()
            print(f"第三个门分解完成，用时 {time4-time3:.4f} 秒, 门序列是 {FT_circuit}")

            # 分解总结
            print("="*50)
            print(f"分解完成，总用时 {time4-time1:.4f} 秒")
            #print("分解完成的量子门是：")
            #print(decompos_gate)
            # 显示矩阵元素
            D = [
                [decompos_gate[0, 0].evalf(), decompos_gate[0, 1].evalf()],
                [decompos_gate[1, 0].evalf(), decompos_gate[1, 1].evalf()],
            ]
            print("\n分解后的量子门矩阵:")
            for row in D:
                print(f"[ {row[0]:.4f}, {row[1]:.4f} ]")

            # 计算门距离
            print("="*50)
            print("门距离:")
            distance = gate_distance_single(D, self.target_matrix)
            print(f"{distance:.4e}")
            print("="*50)
        else:
            phi=parameters["phi"]
            FT_circuit=self.QuantumCompile_1qubit(phi)+"X"
        delta=parameters[0]
        FT_circuit=f"w^{delta}"+FT_circuit
        
        return FT_circuit
    
precise=0.1
target_matrix=UnitaryMatrix([1/sp.sqrt(2),1/sp.sqrt(2),1/sp.sqrt(2),-1/sp.sqrt(2)])
new_compile=NTA(target_matrix,precise)
print("F-T circuit is:".format(new_compile.quantum_complie()))