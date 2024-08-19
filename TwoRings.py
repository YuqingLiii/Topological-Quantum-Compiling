import sympy as sp

class  RingCycIntegers:
    def __init__(self, coeffs):
        """
        初始化一个 Z[ω] 的元素。
        coeffs: 系数列表，表示 a + bω + cω^2 + dω^3
        """
        self.coeffs = coeffs + [0] * (4 - len(coeffs))  # 确保有4项

    def get_str(self):
        terms = []
        for i, coeff in enumerate(self.coeffs):
            if coeff != 0:
                if i == 0:
                    terms.append(f"{coeff}")
                elif i == 1:
                    terms.append(f"{coeff}ω")
                elif i == 2:
                    terms.append(f"{coeff}ω^2")
                elif i == 3:
                    terms.append(f"{coeff}ω^3")
        return " + ".join(terms) if terms else "0"

    def add(self, other):
        result=[0]*4
        for i in range(4):
                result[i]=self.coeffs[i] + other.coeffs[i]

        return RingCycIntegers(result)

    def sub(self, other):
        result=[0]*4
        for i in range(4):
                result[i]=self.coeffs[i] - other.coeffs[i]

        return RingCycIntegers(result)


    def mul(self, other):
        result_coeffs = [0] * 7

        for i, a in enumerate(self.coeffs):
            for j, b in enumerate(other.coeffs):
                result_coeffs[i + j] += a * b
        #print(result_coeffs)
        result_coeffs=self.reduce(result_coeffs)

        return RingCycIntegers(result_coeffs) 
    
    @staticmethod
    def reduce(coeffs):
        """
        将高于 ω^3 的幂次项转换为 ω^3 以下的项。
        """
        if coeffs[4]!=0:
            coeffs[0] += -coeffs[4]  # -1
            coeffs[1] += coeffs[4]   # ω
            coeffs[2] += -coeffs[4]  # -ω^2
            coeffs[3] += coeffs[4]   # ω^3
            coeffs[4] = 0  # 将 ω^4 项归零

        if coeffs[5]!=0:
            coeffs[0] += -coeffs[5]  # -1
            coeffs[5] = 0  # 将 ω^5 项归零
            
        if coeffs[6]!=0:
            coeffs[1]+=-coeffs[6]
            coeffs[6]=0
        return coeffs[:4]
    
    def pow(self,exponent):
        result = RingCycIntegers([1,0,0,0])
        
        base = self
        
        while exponent > 0:
            if exponent % 2 == 1:
                result = result.mul(base)
            base = base.mul(base)
            exponent //= 2
        
        return result
    
    def abs_quadratic(self,is_SubRing=True):
        a=self.get_coeff(0)
        b=self.get_coeff(1)
        c=self.get_coeff(2)
        d=self.get_coeff(3)
        if is_SubRing==True:
            return SubRealRing([a**2+b**2+c**2+d**2+a*b+b*c+c*d,a*b+a*c-a*d+b*c+b*d+c*d])
        else:
            return RingCycIntegers([a**2+b**2+c**2+d**2+a*b+b*c+c*d,0,a*b+a*c-a*d+b*c+b*d+c*d,-(a*b+a*c-a*d+b*c+b*d+c*d)])
        
    def automorphism(self):
        a=self.get_coeff(0)
        b=self.get_coeff(1)
        c=self.get_coeff(2)
        d=self.get_coeff(3)
        return RingCycIntegers([a+d,-c-d,d,b-d])
        
        
    def minus(self):
        a=self.get_coeff(0)
        b=self.get_coeff(1)
        c=self.get_coeff(2)
        d=self.get_coeff(3)
        return RingCycIntegers([a+d,-c-d,d,b-d])
        
    def conjugate_(self):
        a=self.get_coeff(0)
        b=self.get_coeff(1)
        c=self.get_coeff(2)
        d=self.get_coeff(3)
        return RingCycIntegers([-a,-b,-c,-d])
        
    def get_coeff(self, degree):
        """
        返回指定次数的系数。如果次数超出多项式的范围，返回 0。
        """
        if degree < 0 or degree > 3:
            return 0
        return self.coeffs[degree]

    def get_value(self):
        omega=sp.exp(sp.I*sp.pi/5)
        value=self.coeffs[0]+self.coeffs[1]*omega+self.coeffs[2]*omega**2+self.coeffs[3]*omega**3
        return value.evalf()

class SubRealRing:
    def __init__(self, coeffs):
        """
        初始化一个 Z[τ] 的元素。
        coeffs: 系数列表，表示 a + bτ
        """
        self.coeffs = coeffs + [0] * (2 - len(coeffs))  # 确保有2项
        assert self.coeffs[0]%1==0 and self.coeffs[1]%1==0,"环的系数需要是整数"
        self.coeffs[0]=int(self.coeffs[0])
        self.coeffs[1]=int(self.coeffs[1])
    def get_str(self):
        terms = []
        for i, coeff in enumerate(self.coeffs):
            if coeff != 0:
                if i == 0:
                    terms.append(f"{coeff}")
                elif i == 1:
                    terms.append(f"{coeff}τ")
        return " + ".join(terms) if terms else "0"

    def add(self, other):
        result=[0]*2
        for i in range(2):
                result[i]=self.coeffs[i] + other.coeffs[i]
        return SubRealRing(result)

    def sub(self, other):
        result=[0]*2
        for i in range(2):
                result[i]=self.coeffs[i] - other.coeffs[i]
        return SubRealRing(result)

    def transform_ringcycinteger(self):
        return RingCycIntegers([self.coeffs[0],0,self.coeffs[1],-self.coeffs[1]])

    def mul(self, other):
        result_coeffs = [0] * 3

        for i, a in enumerate(self.coeffs):
            for j, b in enumerate(other.coeffs):
                result_coeffs[i + j] += a * b

        result_coeffs=self.reduce(result_coeffs)

        return SubRealRing(result_coeffs) 
    
    @staticmethod
    def reduce(coeffs):
        """
        将高于 τ^1 的幂次项转换为 τ^1 以下的项。
        """
        if coeffs[2]!=0:
            coeffs[0]+=coeffs[2]
            coeffs[1]-=coeffs[2]
            
        return coeffs[:2]
    
    def pow(self,exponent):
        result = SubRealRing([1,0])
        
        base = self
        
        while exponent > 0:
            if exponent % 2 == 1:
                result = result.mul(base)
            base = base.mul(base)
            exponent //= 2
        
        return result       
    
    def automorphism(self):
        result=[0]*2
        result[0]=self.get_coeff(0)-self.get_coeff(1)
        result[1]=-self.get_coeff(1)
        return SubRealRing(result)
    
    def get_coeff(self, degree):
        """
        返回指定次数的系数。如果次数超出多项式的范围，返回 0。
        """
        if degree < 0 or degree > 3:
            return 0
        return self.coeffs[degree]
    
    def N_map(self):
        
        return self.mul(self.automorphism())
    
    def get_value(self):
        
        result = self.coeffs[0] + (self.coeffs[1]) * (sp.sqrt(5)-1)/2

        return result.evalf()
    
    
def test_ring_cyc_integers():
    a = RingCycIntegers([1, 2, 3, 4])
    b = RingCycIntegers([4, 3, 2, 1])

    print("a = ", a.get_str())
    print("b = ", b.get_str())

    # 测试加法
    sum_ab = a.add(b)
    print("a + b = ", sum_ab.get_str())

    # 测试减法
    diff_ab = a.sub(b)
    print("a - b = ", diff_ab.get_str())

    # 测试乘法
    mul_ab = a.mul(b)
    print("a * b = ", mul_ab.get_str())

    # 测试乘方
#    a_pow = a.pow(20)
#    print("a^20 = ", a_pow.get_str())

    # 测试 abs_quadratic
    abs_quad = a.abs_quadratic()
    print("abs_quadratic(a) = ", abs_quad.get_str())

    # 测试 get_value
    a_value = a.get_value()
    print("a 的具体数值 = ", a_value)
  

