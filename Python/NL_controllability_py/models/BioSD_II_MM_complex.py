import sympy as sym
# "BioSD-II" model by Alexis et al, iScience 2021 
# Originally published in: Alexis et al, iScience 2021 
# Corresponds to the model in Fig. 1D

# 3 states
x1 = sym.Symbol('x1')
x2 = sym.Symbol('x2')
x3 = sym.Symbol('x3')
x = [[x1], [x2], [x3]]
# 1 output
h = [x1]

# 1 known input
uu = sym.Symbol('uu')
u = [uu]

# 7 unknown parameters
p1 = sym.Symbol('p1')
p2 = sym.Symbol('p2')
p3 = sym.Symbol('p3')
p4 = sym.Symbol('p4')
p5 = sym.Symbol('p5')
p6 = sym.Symbol('p6')
p7 = sym.Symbol('p7')
p8 = sym.Symbol('p8')
p9 = sym.Symbol('p9')
p = [[p1], [p2], [p3], [p4], [p5], [p6], [p7],[p8],[p9]]

# dynamic equations
f = [[p1*uu+p2-p3*x1*x2-(p4+p5)*x1], [p6*x1*1/(x1+p7)-p8*x2*x3-p5*x2], [p9-p8*x2*x3-p5*x3]]


variables_locales = locals().copy()