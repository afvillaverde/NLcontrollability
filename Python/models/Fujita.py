import sympy as sym

# Model 'Fujita_SciSignal2010' from the benchmarks collection

# 9 states
EGFR = sym.Symbol('EGFR')
pEGFR = sym.Symbol('pEGFR')
pEGFR_Akt = sym.Symbol('pEGFR_Akt')
Akt = sym.Symbol('Akt')
pAkt = sym.Symbol('pAkt')
S6 = sym.Symbol('S6')
pAkt_S6 = sym.Symbol('pAkt_S6')
pS6 = sym.Symbol('pS6')
EGF_EGFR = sym.Symbol('EGF_EGFR')
x = [[EGFR], [pEGFR],[pEGFR_Akt],[Akt],[pAkt],[S6],[pAkt_S6],[pS6],[EGF_EGFR]]


# 1 known input
EGF_u = sym.Symbol('EGF_u')
u = [EGF_u]

# 0 unknown inputs:
w = [];

# 16 unknown parameters
scaleFactor_pEGFR = sym.Symbol('scaleFactor_pEGFR')
scaleFactor_pAkt = sym.Symbol('scaleFactor_pAkt')
scaleFactor_pS6 = sym.Symbol('scaleFactor_pS6')
EGFR_turnover = sym.Symbol('EGFR_turnover')
reaction_1_k1 = sym.Symbol('reaction_1_k1')
reaction_1_k2 = sym.Symbol('reaction_1_k2')
reaction_2_k1 = sym.Symbol('reaction_2_k1')
reaction_2_k2 = sym.Symbol('reaction_2_k2')
reaction_3_k1 = sym.Symbol('reaction_3_k1')
reaction_4_k1 = sym.Symbol('reaction_4_k1')
reaction_5_k1 = sym.Symbol('reaction_5_k1')
reaction_5_k2 = sym.Symbol('reaction_5_k2')
reaction_6_k1 = sym.Symbol('reaction_6_k1')
reaction_7_k1 = sym.Symbol('reaction_7_k1')
reaction_8_k1 = sym.Symbol('reaction_8_k1')
reaction_9_k1 = sym.Symbol('reaction_9_k1')
p = [[scaleFactor_pEGFR], [scaleFactor_pAkt], [scaleFactor_pS6], [EGFR_turnover],
     [reaction_1_k1], [reaction_1_k2], [reaction_2_k1], [reaction_2_k2],
     [reaction_3_k1], [reaction_4_k1], [reaction_5_k1], [reaction_5_k2],
     [reaction_6_k1], [reaction_7_k1], [reaction_8_k1], [reaction_9_k1]]


# dynamic equations
f = [
    [68190.0 * EGFR_turnover - EGF_u * EGFR * reaction_1_k1 + EGF_EGFR * reaction_1_k2 - EGFR * EGFR_turnover],
    [EGF_EGFR * reaction_9_k1 - pEGFR * reaction_4_k1 + pEGFR_Akt * reaction_2_k2 + pEGFR_Akt * reaction_3_k1 - Akt * pEGFR * reaction_2_k1],
    [Akt * pEGFR * reaction_2_k1 - pEGFR_Akt * reaction_3_k1 - pEGFR_Akt * reaction_2_k2],
    [pAkt * reaction_7_k1 + pEGFR_Akt * reaction_2_k2 - Akt * pEGFR * reaction_2_k1],
    [pAkt_S6 * reaction_5_k2 - pAkt * reaction_7_k1 + pAkt_S6 * reaction_6_k1 + pEGFR_Akt * reaction_3_k1 - S6 * pAkt * reaction_5_k1],
    [pAkt_S6 * reaction_5_k2 + pS6 * reaction_8_k1 - S6 * pAkt * reaction_5_k1],
    [S6 * pAkt * reaction_5_k1 - pAkt_S6 * reaction_6_k1 - pAkt_S6 * reaction_5_k2],
    [pAkt_S6 * reaction_6_k1 - pS6 * reaction_8_k1],
    [EGF_u * EGFR * reaction_1_k1 - EGF_EGFR * reaction_1_k2 - EGF_EGFR * reaction_9_k1]
]

# 3 outputs
h = [ [scaleFactor_pEGFR*(pEGFR + pEGFR_Akt)],[
    scaleFactor_pAkt*(pAkt + pAkt_S6)],[
    pS6*scaleFactor_pS6]]
#known_ics = [0,1,1,0,1,0,1,1,1]; 
variables_locales = locals().copy()