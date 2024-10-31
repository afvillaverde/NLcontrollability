% Modelo SIRCD nº 21 paper Gemma COVID
clear;

% 2 states
syms S I R C D
x = [S I R C D];

% 1 output
h = [D];

% one input
syms q;
u = q;

% 4 unknown parameters 
syms beta p1 r mu N
p =[];

% dynamic equations
f = [-beta*S*I/N-q*S+p1*C; 
    beta*S*I/N-(r+mu)*I; 
    r*I; 
    q*S-p1*C;
    mu*I];

% initial conditions
ics  = []; 
known_ics = [0,0,0,0,0];

save('SIRCD_knownparam','x','p','h','f','u','ics','known_ics');

p=[beta,p1,r,mu];

save('SIRCD_unknownparam','x','p','h','f','u','ics','known_ics');

p=[beta,p1,r];

save('SIRCD_unknownparam_mu','x','p','h','f','u','ics','known_ics');


f=[-beta*S*I/N-q*S+p1*C;
    beta*S*I/N-(r+1)*I;
    r*I;
    q*S-p1*C;
    I
    ];

save('SIRCD_unknownparam_repar','x','p','h','f','u','ics','known_ics');