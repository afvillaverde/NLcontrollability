from sympy import symbols,Matrix
ics = []
known_ics = []
w = []
x1, x2, x3, x4, x5 = symbols('x1, x2, x3, x4, x5')
x = [[x1], [x2], [x3], [x4], [x5]]
Kg, Ka, l, rgupP, raupP, raoverP, kpts, Ka, kacs, rgupC, raoverC, raupC, D, yh, kdeg = symbols('Kg, Ka, l, rgupP, raupP, raoverP, kpts, Ka, kacs, rgupC, raoverC, raupC, D, yh, kdeg')
yg, ya, kover, kg, ka, n, m, Oa, Og = symbols('yg, ya, kover, kg, ka, n, m, Oa, Og')
p = [[yg], [ya], [kover], [kg], [ka], [n], [m], [Oa], [Og]]
h = [[x1], [x2], [x3 + x5]]
Gin = symbols('Gin')
u = [Gin]
f = [[D*(Gin - x1) - Oa**n*kg*x1*x3/((Kg + x1)*(Oa**n + x2**n)) - Oa**n*kpts*x1*x4/((Kg + x1)*(Oa**n + x2**n))], [-D*x2 + x3*(-Og**m*ka*x2/((Ka + x2)*(Og**m + (Oa**n*kg*x1/((Kg + x1)*(Oa**n + x2**n)))**m)) + kover*(Oa**n*kg*x1/((Kg + x1)*(Oa**n + x2**n)) - l)) + x4*(-Og**m*ka*x2/((Ka + x2)*(Og**m + (Oa**n*kpts*x1/((Kg + x1)*(Oa**n + x2**n)))**m)) - kacs*x2/(kacs + x2) + kover*(Oa**n*kpts*x1/((Kg + x1)*(Oa**n + x2**n)) - l))], [-D*x3 - kdeg*x3 + x3*(1 - yh)*(Oa**n*kg*x1*yg/((Kg + x1)*(Oa**n + x2**n)) + ya*(Og**m*ka*x2/((Ka + x2)*(Og**m + (Oa**n*kg*x1/((Kg + x1)*(Oa**n + x2**n)))**m)) - kover*(Oa**n*kg*x1/((Kg + x1)*(Oa**n + x2**n)) - l)))], [-D*x4 - kdeg*x4 + x4*(Oa**n*kpts*x1*yg/((Kg + x1)*(Oa**n + x2**n)) + ya*(Og**m*ka*x2/((Ka + x2)*(Og**m + (Oa**n*kpts*x1/((Kg + x1)*(Oa**n + x2**n)))**m)) + kacs*x2/(kacs + x2) - kover*(Oa**n*kpts*x1/((Kg + x1)*(Oa**n + x2**n)) - l)))], [-D*x5 - kdeg*x5 + x3*yh*(Oa**n*kg*x1*yg/((Kg + x1)*(Oa**n + x2**n)) + ya*(Og**m*ka*x2/((Ka + x2)*(Og**m + (Oa**n*kg*x1/((Kg + x1)*(Oa**n + x2**n)))**m)) - kover*(Oa**n*kg*x1/((Kg + x1)*(Oa**n + x2**n)) - l)))]]
fu = Matrix([[D*(Gin - x1)/Gin], [0], [0], [0], [0]])
hu = Matrix([[0], [0], [0]])
fxw = Matrix([[Oa**n*x1*(-kg*x3 - kpts*x4)/((Kg + x1)*(Oa**n + x2**n))], [-D*x2 - x3*(Og**m*ka*x2/((Ka + x2)*(Og**m + (Oa**n*kg*x1/((Kg + x1)*(Oa**n + x2**n)))**m)) - kover*(Oa**n*kg*x1/((Kg + x1)*(Oa**n + x2**n)) - l)) - x4*(Og**m*ka*x2/((Ka + x2)*(Og**m + (Oa**n*kpts*x1/((Kg + x1)*(Oa**n + x2**n)))**m)) + kacs*x2/(kacs + x2) - kover*(Oa**n*kpts*x1/((Kg + x1)*(Oa**n + x2**n)) - l))], [-D*x3 - kdeg*x3 - x3*(yh - 1)*(Oa**n*kg*x1*yg/((Kg + x1)*(Oa**n + x2**n)) + ya*(Og**m*ka*x2/((Ka + x2)*(Og**m + (Oa**n*kg*x1/((Kg + x1)*(Oa**n + x2**n)))**m)) - kover*(Oa**n*kg*x1/((Kg + x1)*(Oa**n + x2**n)) - l)))], [-D*x4 - kdeg*x4 + x4*(Oa**n*kpts*x1*yg/((Kg + x1)*(Oa**n + x2**n)) + ya*(Og**m*ka*x2/((Ka + x2)*(Og**m + (Oa**n*kpts*x1/((Kg + x1)*(Oa**n + x2**n)))**m)) + kacs*x2/(kacs + x2) - kover*(Oa**n*kpts*x1/((Kg + x1)*(Oa**n + x2**n)) - l)))], [-D*x5 - kdeg*x5 + x3*yh*(Oa**n*kg*x1*yg/((Kg + x1)*(Oa**n + x2**n)) + ya*(Og**m*ka*x2/((Ka + x2)*(Og**m + (Oa**n*kg*x1/((Kg + x1)*(Oa**n + x2**n)))**m)) - kover*(Oa**n*kg*x1/((Kg + x1)*(Oa**n + x2**n)) - l)))]])
hxw = Matrix([[x1], [x2], [x3 + x5]])
variables_locales = locals().copy()