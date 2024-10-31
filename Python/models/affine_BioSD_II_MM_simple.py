from sympy import symbols,Matrix
ics = []
known_ics = []
w = []
x1, x2, x3 = symbols('x1, x2, x3')
x = [[x1], [x2], [x3]]
p1, p2, p3, p4, p5, p6, p7, p8 = symbols('p1, p2, p3, p4, p5, p6, p7, p8')
p = [[p1], [p2], [p3], [p4], [p5], [p6], [p7], [p8]]
h = [x1]
uu = symbols('uu')
u = [uu]
f = [[p1*uu + p2 - p3*x1*x2 - p4*x1], [p5*x1/(p6 + x1) - p7*x2*x3], [-p7*x2*x3 + p8]]
fu = Matrix([[p1], [0], [0]])
hu = Matrix([[0]])
fxw = Matrix([[p2 - p3*x1*x2 - p4*x1], [(p5*x1 - p7*x2*x3*(p6 + x1))/(p6 + x1)], [-p7*x2*x3 + p8]])
hxw = Matrix([[x1]])
variables_locales = locals().copy()