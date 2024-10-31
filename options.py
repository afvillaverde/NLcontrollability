#==============================================================================
#THE USER CAN DEFINE THE PROBLEM AND SET OPTIONS IN THE FOLLOWING LINES:
#==============================================================================
from math import inf
###############################################################################
# (1) NAME OF THE MODEL TO BE STUDIED:
modelname = '1D_BIG' #name of a .py file placed in 'models' folder
# (2) ANALYSIS OPTIONS:
LCtest = 0 # 1 if you want to check the Linearization Condition
              #0 if not

numericLC=0   # 0 if you want to check the Linearization Condition symbolicaly
              #1 for numeric computation
              
ARC=0        # 1 if you want to check the Accesibility Rank Condition
              #0 if not
              
LARCtest=1    # 1 if you want to check the Lie Algebraic Rank Condition
              #0 if not

GSC=0        # 1 if you want to check the Sussmanns General Sufficient Condition
              #0 if not

maxtime = 5000  # max. time for each test (seconds)

x0= [1,1,1] #Specific initial point. If no  point is given (x0="-") it 
                           #would try to compute equilibrium points.

