from sympy import symbols,Matrix
ics = []
known_ics = []
w = []
EGFR, pEGFR, pEGFR_Akt, Akt, pAkt, S6, pAkt_S6, pS6, EGF_EGFR = symbols('EGFR, pEGFR, pEGFR_Akt, Akt, pAkt, S6, pAkt_S6, pS6, EGF_EGFR')
x = [[EGFR], [pEGFR], [pEGFR_Akt], [Akt], [pAkt], [S6], [pAkt_S6], [pS6], [EGF_EGFR]]
scaleFactor_pEGFR, scaleFactor_pAkt, scaleFactor_pS6, EGFR_turnover, reaction_1_k1, reaction_1_k2, reaction_2_k1, reaction_2_k2, reaction_3_k1, reaction_4_k1, reaction_5_k1, reaction_5_k2, reaction_6_k1, reaction_7_k1, reaction_8_k1, reaction_9_k1 = symbols('scaleFactor_pEGFR, scaleFactor_pAkt, scaleFactor_pS6, EGFR_turnover, reaction_1_k1, reaction_1_k2, reaction_2_k1, reaction_2_k2, reaction_3_k1, reaction_4_k1, reaction_5_k1, reaction_5_k2, reaction_6_k1, reaction_7_k1, reaction_8_k1, reaction_9_k1')
p = [[scaleFactor_pEGFR], [scaleFactor_pAkt], [scaleFactor_pS6], [EGFR_turnover], [reaction_1_k1], [reaction_1_k2], [reaction_2_k1], [reaction_2_k2], [reaction_3_k1], [reaction_4_k1], [reaction_5_k1], [reaction_5_k2], [reaction_6_k1], [reaction_7_k1], [reaction_8_k1], [reaction_9_k1]]
h = [[scaleFactor_pEGFR*(pEGFR + pEGFR_Akt)], [scaleFactor_pAkt*(pAkt + pAkt_S6)], [pS6*scaleFactor_pS6]]
EGF_u = symbols('EGF_u')
u = [EGF_u]
f = [[-EGFR*EGFR_turnover - EGFR*EGF_u*reaction_1_k1 + 68190.0*EGFR_turnover + EGF_EGFR*reaction_1_k2], [-Akt*pEGFR*reaction_2_k1 + EGF_EGFR*reaction_9_k1 - pEGFR*reaction_4_k1 + pEGFR_Akt*reaction_2_k2 + pEGFR_Akt*reaction_3_k1], [Akt*pEGFR*reaction_2_k1 - pEGFR_Akt*reaction_2_k2 - pEGFR_Akt*reaction_3_k1], [-Akt*pEGFR*reaction_2_k1 + pAkt*reaction_7_k1 + pEGFR_Akt*reaction_2_k2], [-S6*pAkt*reaction_5_k1 - pAkt*reaction_7_k1 + pAkt_S6*reaction_5_k2 + pAkt_S6*reaction_6_k1 + pEGFR_Akt*reaction_3_k1], [-S6*pAkt*reaction_5_k1 + pAkt_S6*reaction_5_k2 + pS6*reaction_8_k1], [S6*pAkt*reaction_5_k1 - pAkt_S6*reaction_5_k2 - pAkt_S6*reaction_6_k1], [pAkt_S6*reaction_6_k1 - pS6*reaction_8_k1], [EGFR*EGF_u*reaction_1_k1 - EGF_EGFR*reaction_1_k2 - EGF_EGFR*reaction_9_k1]]
fu = Matrix([[-EGFR*reaction_1_k1], [0], [0], [0], [0], [0], [0], [0], [EGFR*reaction_1_k1]])
hu = Matrix([[0], [0], [0]])
fxw = Matrix([[-EGFR*EGFR_turnover + 68190.0*EGFR_turnover + EGF_EGFR*reaction_1_k2], [-Akt*pEGFR*reaction_2_k1 + EGF_EGFR*reaction_9_k1 - pEGFR*reaction_4_k1 + pEGFR_Akt*reaction_2_k2 + pEGFR_Akt*reaction_3_k1], [Akt*pEGFR*reaction_2_k1 - pEGFR_Akt*reaction_2_k2 - pEGFR_Akt*reaction_3_k1], [-Akt*pEGFR*reaction_2_k1 + pAkt*reaction_7_k1 + pEGFR_Akt*reaction_2_k2], [-S6*pAkt*reaction_5_k1 - pAkt*reaction_7_k1 + pAkt_S6*reaction_5_k2 + pAkt_S6*reaction_6_k1 + pEGFR_Akt*reaction_3_k1], [-S6*pAkt*reaction_5_k1 + pAkt_S6*reaction_5_k2 + pS6*reaction_8_k1], [S6*pAkt*reaction_5_k1 - pAkt_S6*reaction_5_k2 - pAkt_S6*reaction_6_k1], [pAkt_S6*reaction_6_k1 - pS6*reaction_8_k1], [EGF_EGFR*(-reaction_1_k2 - reaction_9_k1)]])
hxw = Matrix([[scaleFactor_pEGFR*(pEGFR + pEGFR_Akt)], [scaleFactor_pAkt*(pAkt + pAkt_S6)], [pS6*scaleFactor_pS6]])
variables_locales = locals().copy()