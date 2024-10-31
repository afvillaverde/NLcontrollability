import sympy as sym
# "BioSD-I" model by Alexis et al, iScience 2021 
# Originally published in: Alexis et al, iScience 2021 
# Corresponds to the model in Fig. 1C

# 2 states
x1 = sym.Symbol('x1')
x2 = sym.Symbol('x2')
x = [[x1], [x2]]
# 1 output
h = [x1]

# 1 known input
uu = sym.Symbol('uu')
u = [uu]

# 6 unknown parameters
p1 = sym.Symbol('p1')
p2 = sym.Symbol('p2')
p3 = sym.Symbol('p3')
p4 = sym.Symbol('p4')
p5 = sym.Symbol('p5')
p6 = sym.Symbol('p6')
p = [[p1], [p2], [p3], [p4], [p5], [p6]]

# dynamic equations
f = [[p1*uu+p2-p3*x1*x2-p4*x2], [p6*x1-p5]]


variables_locales = locals().copy()