import sympy as sym
#ctrl_ex_61
#4 states
x1 = sym.Symbol('x1')
x2 = sym.Symbol('x2')
x3 = sym.Symbol('x3')
x4 = sym.Symbol('x4')
x = [[x1], [x2], [x3], [x4]]

# 2 output

h = [[x1],[x2]]

# 1 known input
uu= sym.Symbol('uu')
u = [uu]

# 0 unknown parameters
p = []

# dynamic equations
f = [[uu], [x1], [x1**3], [-x2**7+x3**2]]


variables_locales = locals().copy()