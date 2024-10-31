% Models on paper Control-theoretic analysis of synthetic biology circuits with differential geometry
clear;

% 2 states
syms x1 x2
x = [x1; x2];

% 1 input
syms u;

% 6 parameters 
syms p1 p2 p3 p4 p5 p6 
p =[p1; p2; p3; p4; p5; p6];

% initial conditions
ics  = []; 
known_ics = [0,0];

% dynamic equations
f = [p1*u+p2-p3*x1*x2-p4*x1;
    p6*x1-p5];

% 1 output
h = x1 ;

save('BioSD_I','x','p','u','h','f','ics','known_ics');

clear;

% 3 states
syms x1 x2 x3
x = [x1; x2; x3];

% 1 input
syms u;

% 7 parameters 
syms p1 p2 p3 p4 p5 p6 p7
p =[p1; p2; p3; p4; p5; p6; p7];

% initial conditions
ics  = []; 
known_ics = [0,0];

% dynamic equations
f = [p1*u+p2-p3*x1*x2-p4*x1;
    p6*x1-p7*x2*x3;
    p5-p7*x2*x3];

% 1 output
h = x1 ;

save('BioSD_II','x','p','u','h','f','ics','known_ics');

clear;

% 3 states
syms x1 x2 x3
x = [x1; x2; x3];

% 1 input
syms u;

% 7 parameters 
syms p1 p2 p3 p4 p5 p6 p7
p =[p1; p2; p3; p4; p5; p6; p7];

% initial conditions
ics  = []; 
known_ics = [0,0];

% dynamic equations
f = [p1*u+p2-p3*x1*x2-p4*x1+p3*x1*x3;
    p6*x1-p7*x2*x3;
    p5-p7*x2*x3];

% 1 output
h = x1 ;

save('BioSD_III','x','p','u','h','f','ics','known_ics');

clear;

% 3 states
syms x1 x2 x3
x = [x1; x2; x3];

% 1 input
syms u;

% 8 parameters 
syms p1 p2 p3 p4 p5 p6 p7 p8
p =[p1; p2; p3; p4; p5; p6; p7; p8];

% initial conditions
ics  = []; 
known_ics = [0,0];

% dynamic equations
f = [p1*u+p2-p3*x1*x2-p4*x1;
    p5*x1/(x1+p6)-p7*x2*x3;
    p8-p7*x2*x3];

% 1 output
h = x1 ;

save('BioSD_II_MM_simple','x','p','u','h','f','ics','known_ics');

clear;

% 3 states
syms x1 x2 x3
x = [x1; x2; x3];

% 1 input
syms u;

% 8 parameters 
syms p1 p2 p3 p4 p5 p6 p7 p8 p9
p =[p1; p2; p3; p4; p5; p6; p7; p8; p9];

% initial conditions
ics  = []; 
known_ics = [0,0];

% dynamic equations
f = [p1*u+p2-p3*x1*x2-(p4+p5)*x1;
    p6*x1/(x1+p7)-p8*x2*x3-p5*x2;
    p9-p8*x2*x3-p5*x3];

% 1 output
h = x1 ;

save('BioSD_II_MM_complex','x','p','u','h','f','ics','known_ics');

clear;

% 3 states
syms x1 x2 x3
x = [x1; x2; x3];

% 8 parameters 
syms p1 p2 p3 p4 p5 p6 p7 p8 p9
p =[ p2; p3; p4; p5; p6; p7; p8; p9];

% 1 input
u=p1;

% initial conditions
ics  = []; 
known_ics = [0,0];

% dynamic equations
f = [p1-p2*x1-p3*x1+p4*(p1/p2-x1)*x2+p5*(p1/p2-x1)*x3;
    p6-p2*x2-p4*(p1/p2-x1)*x2+p7*x1*(p6/p2-x2);
    p8-p2*x3-p5*(p1/p2-x1)*x3+p9*x1*(p8/p2-x3)];

% 1 output
h = x1 ;

save('Dichotomous_Feedback_BettaHK','x','p','u','h','f','ics','known_ics');

% 8 parameters 
syms p1 p2 p3 p4 p5 p6 p7 p8 p9
p =[ p1; p2; p3; p4; p5; p7; p8; p9];

% 1 input
u=p6;

save('Dichotomous_Feedback_BettaRR','x','p','u','h','f','ics','known_ics');

% 8 parameters 
syms p1 p2 p3 p4 p5 p6 p7 p8 p9
p =[ p1; p2; p3; p4; p5; p6; p7; p9];

u=p8;

save('Dichotomous_Feedback_BettaSR','x','p','u','h','f','ics','known_ics');

% 8 parameters 
syms p1 p2 p3 p4 p5 p6 p7 p8 p9
p =[ p1; p2; p4; p5; p6; p7; p8; p9];

% 1 input
u=p3;

save('Dichotomous_Feedback_kap','x','p','u','h','f','ics','known_ics');

clear;

% 3 states
syms x1 x2 x3
x = [x1; x2; x3];

% 8 parameters 
syms p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11
p =[p2; p3; p4; p5; p6; p7; p8; p9; p10; p11];

% 1 input
u = p1;

% initial conditions
ics  = []; 
known_ics = [0,0];

% dynamic equations
f = [p1-p2*x1-p10*p3/(p3+p11)*x1+p4*(p1/p2-x1)*x2+p5*(p1/p2-x1)*x3;
    p6-p2*x2-p4*(p1/p2-x1)*x2+p7*x1*(p6/p2-x2);
    p8-p2*x3-p5*(p1/p2-x1)*x3+p9*x1*(p8/p2-x3)];

% 1 output
h = x1 ;

save('Dichotomous_Feedback_I_BettaHK','x','p','u','h','f','ics','known_ics');

% 8 parameters 
syms p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11
p =[p1; p2; p3; p4; p5; p7; p8; p9; p10; p11];

% 1 input
u = p6;

save('Dichotomous_Feedback_I_BettaRR','x','p','u','h','f','ics','known_ics');

% 8 parameters 
syms p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11
p =[p1; p2; p3; p4; p5; p6; p7; p9; p10; p11];

% 1 input
u = p8;

save('Dichotomous_Feedback_I_BettaSR','x','p','u','h','f','ics','known_ics');

% 8 parameters 
syms p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11
p =[p1; p2; p4; p5; p6; p7; p8; p9; p10; p11];

% 1 input
u = p3;

save('Dichotomous_Feedback_I','x','p','u','h','f','ics','known_ics');

clear;

% 4 states
syms t s c T
x = [t; s; c; T];

% 2 input
syms u1 u2;
u = [u1;u2];
% 14 parameters 
syms p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11 p12 p13 p14
p =[p1; p2; p3; p4; p5; p6; p7; p8; p9; p10; p11;p12;p13; p14];

gammaT=(p1)/(1+(T/(p2*(u1/p3)^p4))^p5)+(p14*(T/(p2*((1+(u1/p3))^p4)^p5)))/(1+(T/(p2*((1+(u1/p3))^p4)^p5)));
gammaR=p8*u2/(p9+u2);
% 1 output
h = t ;

% initial conditions
ics  = []; 
known_ics = [0,0,0,0];

% dynamic equations
f = [gammaT-p6*t-p7*t*s;
    gammaR-p10*s-p7*t*s;
    p7*t*s-p11*c;
    p12*t-p13*T];

save('Kelly_1','x','p','u','h','f','ics','known_ics');

clear;

% 4 states
syms t s c T
x = [t; s; c; T];

% 2 input
syms u1 gammaR;

u = [gammaR];
% 14 parameters 
syms p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11 p12 p13 p14 
p =[p1; p2; p3; p4; p5; p6; p7; p10; p11;p12;p13; p14];

gammaT=(p1)/(1+(T/(p2*(u1/p3)^p4))^p5)+(p14*(T/(p2*((1+(u1/p3))^p4)^p5)))/(1+(T/(p2*((1+(u1/p3))^p4)^p5)));
% 1 output
h = t ;

% initial conditions
ics  = []; 
known_ics = [0,0,0,0];

% dynamic equations
f = [gammaT-p6*t-p7*t*s;
    gammaR-p10*s-p7*t*s;
    p7*t*s-p11*c;
    p12*t-p13*T];


save('Kelly_1_gr','x','p','u','h','f','ics','known_ics');

clear;

% 4 states
syms s r c R
x = [r; s; c; R];

% 2 input
syms u2 u3;
u = [u2;u3];

% 11 parameters 
syms p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11
p =[p1; p2; p3; p4; p5; p6; p7; p8; p9; p10; p11];

gammaX_star=(p9*u3)/(u3+p10);
gammaR_star=(p1*R*u2)/((u2+p2)*(p3+(R*u2/u2+p2)));

% initial conditions
ics  = []; 
known_ics = [0,0];

% dynamic equations
f = [gammaX_star-p4*r-p5*r*s;
    gammaR_star-p11*s-p5*r*s;
    p5*r*s-p6*c;
    p7*r-p8*R];

% 1 output
h = [R] ;

save('Kelly_2','x','p','u','h','f','ics','known_ics');

clear;

% 4 states
syms s r c R
x = [r; s; c; R];

% 2 input
syms gammaX_star gammaR_star;
u = [gammaR_star;gammaR_star];

% 11 parameters 
syms p4 p5 p6 p7 p8 p10 p11
p =[ p4; p5; p6; p7; p8; p11];

% initial conditions
ics  = []; 
known_ics = [0,0];

% dynamic equations
f = [gammaX_star-p4*r-p5*r*s;
    gammaR_star-p11*s-p5*r*s;
    p5*r*s-p6*c;
    p7*r-p8*R];

% 1 output
h = [R] ;

save('Kelly_2_gxgr','x','p','u','h','f','ics','known_ics');