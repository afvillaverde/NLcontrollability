import sympy as sym

# Model Microbial Community

# 5 states
x1 = sym.Symbol('x1')
x2 = sym.Symbol('x2')
x3 = sym.Symbol('x3')
x4 = sym.Symbol('x4')
x5 = sym.Symbol('x5')
x = [[x1], [x2],[x3],[x4],[x5]]


# 1 known input
Gin = sym.Symbol('Gin')
u = [Gin]

# 0 unknown inputs:
w = [];

# 9 unknown parameters
yg = sym.Symbol('yg')
ya = sym.Symbol('ya')
kover = sym.Symbol('kover')
kg = sym.Symbol('kg')
ka = sym.Symbol('ka')
n = sym.Symbol('n')
m = sym.Symbol('m')
Oa = sym.Symbol('Oa')
Og = sym.Symbol('Og')
p = [[yg], [ya], [kover], [kg], [ka], [n], [m], [Oa], [Og]]


#Define the following variables as known variables 
#productor
Kg = sym.Symbol('Kg')
Ka = sym.Symbol('Ka')
l = sym.Symbol('l')
rgupP = kg*(x1/(x1+Kg))*((Oa**n)/((x2**n)+(Oa**n)))  
raupP = ka*(x2/(x2+Ka))*((Og**m)/((rgupP**m)+(Og**m)))  
raoverP = kover*(rgupP-l)

#cleaner
kpts = sym.Symbol('kpts')  
Ka = sym.Symbol('Ka')  
kacs = sym.Symbol('kacs')  
rgupC = kpts*(x1/(x1+Kg))*((Oa**n)/((x2**n)+(Oa**n)))  # Corregido
raoverC = kover*(rgupC-l)
raupC = (ka*(x2/(x2+Ka))*((Og**m)/((rgupC**m)+(Og**m))))+(kacs*(x2/(x2+kacs)))
D = sym.Symbol('D')
yh = sym.Symbol('yh')
kdeg= sym.Symbol('kdeg')

# dynamic equations
f =  [
        [-(rgupP)*x3 - rgupC*x4 + D*(Gin - x1)],
        [(raoverP - raupP)*x3 + (raoverC - raupC)*x4 - D*x2],
        [(1 - yh)*(yg*rgupP + ya*(raupP - raoverP))*x3 - kdeg*x3 - D*x3],
        [(yg*rgupC + ya*(raupC - raoverC))*x4 - kdeg*x4 - D*x4],
        [yh*(yg*rgupP + ya*(raupP - raoverP))*x3 - kdeg*x5 - D*x5]
    ]

# 3 outputs
h    =  [[x1],[x2],[x3+x5]];
#known_ics = [1;1;1;1;1]; 
variables_locales = locals().copy()