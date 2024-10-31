from sympy import symbols,Matrix
ics = []
known_ics = []
w = []
x1, x2, x3 = symbols('x1, x2, x3')
x = [[x1], [x2], [x3]]
p1, p2, p3, p4, p5, p6, p7, p8, p9, p11, p12 = symbols('p1, p2, p3, p4, p5, p6, p7, p8, p9, p11, p12')
p = [[p1], [p2], [p3], [p4], [p5], [p6], [p7], [p8], [p9], [p11], [p12]]
h = [x2]
uu = symbols('uu')
u = [uu]
f = [[p1 - p2*x1 - p3*x1 + p4*x2*(p5 - x1) + p6*x3*(p5 - x1)], [-p2*x2 - p4*x2*(p5 - x1) + p7 + p8*x1*(p9 - x2)], [p11*x1*(p12 - x3) - p2*x3 - p6*x3*(p5 - x1) + uu]]
fu = Matrix([[0], [0], [1]])
hu = Matrix([[0]])
fxw = Matrix([[p1 - p2*x1 - p3*x1 + p4*x2*(p5 - x1) + p6*x3*(p5 - x1)], [-p2*x2 - p4*x2*(p5 - x1) + p7 + p8*x1*(p9 - x2)], [p11*x1*(p12 - x3) - p2*x3 - p6*x3*(p5 - x1)]])
hxw = Matrix([[x2]])
variables_locales = locals().copy()