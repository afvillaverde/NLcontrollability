import sympy as sym
# "Dichotomous Feedback" model by Sootla et al, royalsocietypublishing.org/journal/rsif 
# Originally published in: Sootla et al, royalsocietypublishing.org/journal/rsif 
# Corresponds to the model in equations 2.7


#nomenclatura
# HK=x1
# RR=x2
# SR=x3
# Bhk=p1
# delta=p2
# Kap=p3
# kt=p4
# HKtot=p5
# Ktc=p6
# Brr=p7
# Kp=p8
# RRtot=p9
# Bsr=u
# Kpc=p11
# SRtot=p12

# 3 states
x1 = sym.Symbol('x1')
x2 = sym.Symbol('x2')
x3 = sym.Symbol('x3')
x = [[x1], [x2],[x3]]
# 1 output
h = [x2]

# 1 known input
uu = sym.Symbol('uu')
u = [uu]
# 12 unknown parameters
p1 = sym.Symbol('p1')
p2 = sym.Symbol('p2')
p3 = sym.Symbol('p3')
p4 = sym.Symbol('p4')
p5 = sym.Symbol('p5')
p6 = sym.Symbol('p6')
p7 = sym.Symbol('p7')
p8 = sym.Symbol('p8')
p9 = sym.Symbol('p9')
#p10 = sym.Symbol('p10')
p11 = sym.Symbol('p11')
p12 = sym.Symbol('p12')

p = [[p1], [p2], [p3], [p4], [p5], [p6],[p7],[p8],[p9],[p11],[p12]]

# dynamic equations
f = [[p1-p2*x1-p3*x1+p4*(p5-x1)*x2+p6*x3*(p5-x1)], [p7-p2*x2-p4*x2*(p5-x1)+p8*x1*(p9-x2)],[uu-p2*x3-p6*x3*(p5-x1)+p11*x1*(p12-x3)]]


variables_locales = locals().copy()