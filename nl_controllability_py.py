"""
=================================================================================================================
 -> NL controllability 
-----------------------------------------------------------------------------------------------------------------
  * Based on the Matlab NL controllability version by Alejandro Fernandez Villaverde (afvillaverde@uvigo.gal)
-----------------------------------------------------------------------------------------------------------------
  * User input can be specified in file 'options.py' 
=================================================================================================================
"""
def nl_controllability (*args):
    import os, math, random
    import numpy as np
    from sympy import eye,zeros, simplify, nextprime,collect, Matrix,transpose,solve,Eq
    import sympy as sym
    from time import time
    from math import ceil
    from pathlib import Path
    import importlib
    from sympy.matrices.common import ShapeError
    # NL_controllability directory...
    dir = os.path.dirname(os.path.abspath(__file__))
    ########################################################################
    # Create 'results' folder in NL directory in case it does not exist
    path = Path('{}/results'.format(dir))
    path.mkdir(parents=True, exist_ok=True)
    ########################################################################
    print('\n\n ---------------------------------------')
    print(' >>> NL Controllability toolbox for phyton <<<')
    print(' --------------------------------------- \n')
    t_start = time()
    ########################################################################
    if len(args)==0:
        import options
        model = importlib.import_module('models.{}'.format(options.modelname))
        print(' Analyzing the {} model... \n'.format(options.modelname))
    if len(args)==1:
        model = importlib.import_module('models.{}'.format(options.modelname))
        print(' Analyzing the {} model... \n'.format(options.modelname))
    ###################################
    ###################################
    #Check or (convert) the model to the affine formulation:
    try:
        prueba=model.u
    except AttributeError:
        print(f"\n{options.modelname} is an uncontrolled system.")
        print("Its controllability cannot be analysed.\n")
        return
    #############Output Variable initialization
    is_affine=0
    fxw=[]
    fu=[]
    hxw=[]
    hu=[]
    # nº of states
    ns=len(model.x)
    #nº of outputs
    m=len(model.h)
    #nº of known inputs
    nu=len(model.u)
    #nº of unknown inputs
    if hasattr(model,'w'):
        nw=len(model.w)
    #nº of unknown parameters
    np=len(model.p)

    #################################################
    #Checking if the system is affine in known inputs and storing coefficients
    #if the file 'affine' already exists
    file_path = os.path.join(os.getcwd(), 'models',f'affine_{options.modelname}.py')
    if os.path.exists(file_path):
        print(f'\n >>>affine_{options.modelname} file already exists...')
        affine_model = importlib.import_module('models.affine_{}'.format(options.modelname))
        is_affine=1
    #if the file 'affine' doesn't exist
    else:
        print('\n >>> Converting model into affine-in-inputs form...')
        #Initialization of the coefficients of known inputs:
        #Dimension ns*nu (number of states * number of known inputs)
        fu=zeros(ns,nu)
        #Initialization of the coefficients of known outputs:
        #Dimension m*nu (number of outputs* nº number of known inputs)
        hu=zeros(m,nu)
        for j in range(0,nu):
            uj=(model.u)[j]
            #states
            for i in range(0,ns):
                expression = (model.f)[i][0]
                print(expression)
                terms= expression.as_ordered_terms()
                result=0
                for term in terms:
                    if term.has(uj):
                        result+=term
                simplified_expression = simplify(result)
                if simplified_expression==uj:
                    fu[i]=1
                else:                             
                    if simplified_expression!=0:
                        factored_expr = collect(simplified_expression, uj)
                        rearranged_expr = factored_expr/ uj
                        fu[i]=rearranged_expr
            #outputs
            for i in range (0,m):
                if hasattr(model,'w') and (len(model.w)>0):
                    hu[i]=model.w[i]
                else:
                    hu[i]=0
        # State-unknown inputs 0-augmented system dynamics:            
        u=Matrix(model.u)
        f=Matrix(model.f)
        h=Matrix(model.h)
        if isinstance(u, Matrix) and u.cols > 1:
            u=transpose(u)
        try:
            fxw=simplify(f-fu*u)
        except ShapeError:
            print(f'\n {options.modelname} is an uncontrolled system. ');
            print('Its controllability cannot be analysed. \n\n')
            return
       # State-unknown inputs output dynamics:
        hxw=simplify(h-hu*u)
       #######Save affine model
        affine_modelname = f'affine_{options.modelname}'
        fullaffinename = os.path.join(os.getcwd(), 'models', affine_modelname)
        phyton_code = "from sympy import symbols,Matrix\n"
        # Check variable existence
        ics = [] if 'ics' not in locals() else ics
        known_ics = [] if 'known_ics' not in locals() else known_ics
        w = [] if 'w' not in locals() else w
        phyton_code += f"ics = {ics}\n"
        phyton_code += f"known_ics = {known_ics}\n"
        phyton_code += f"w = {w}\n"   
        #other variables_to_assign=['x','u','p', 'f', 'h']
        x=model.x
        u=model.u
        p=model.p
        f=model.f
        h=model.h
        #saving x
        print(model.x)
        length_x = len(x)
        x_new = [f"{x[i]}" for i in range(length_x)]
        x_new_str = ', '.join(x_new)
        x_new_str = x_new_str.replace("[", "").replace("]", "")
        phyton_code += f"{x_new_str} = symbols('{x_new_str}')\n"
        x_completo = [sym.Symbol(f'{x[i]}') for i in range(len(x))]
        phyton_code += f"x = {x_completo}\n"      
        #Saving known variables
        # Ask the user if the model has known variables
        tiene_variables_conocidas = input("Does the model have known variables? (yes/no): ").lower()
        if tiene_variables_conocidas == "yes":
            #ask the user for a vector with the known variables
            vector_variables_conocidas = input("Enter a vector with the known variables separated by commas: ")
            variables_conocidas = vector_variables_conocidas.split(',')
            variablesK_str = ', '.join(variables_conocidas)
            phyton_code += f"{variablesK_str} = symbols('{variablesK_str}')\n"     
        else:
            print("The model has no known variables.")        
        #saving p
        length_p = len(p)
        if length_p>0:
            p_new = [f"{p[i]}" for i in range(length_p)]
            p_new_str = ', '.join(p_new)
            p_new_str = p_new_str.replace("[", "").replace("]", "")
            phyton_code += f"{p_new_str} = symbols('{p_new_str}')\n"
            p_completo = [sym.Symbol(f'{p[i]}') for i in range(len(p))]
            phyton_code += f"p = {p_completo}\n"
        #saving h
        length_h = len(h)
        h_completo = [sym.Symbol(f'{h[i]}') for i in range(len(h))]
        phyton_code += f"h = {h_completo}\n"
        #saving u0, model 1D_BIG
        try:
            print (model.u0)
            phyton_code += "u0 = symbols('u0')\n"
        except:
            1-1
        #saving u
        length_u = len(u)
        u_new = [f"{u[i]}" for i in range(length_u)]
        u_new_str = ', '.join(u_new)
        u_new_str = u_new_str.replace("[", "").replace("]", "")
        phyton_code += f"{u_new_str} = symbols('{u_new_str}')\n"
        u_completo = [sym.Symbol(f'{u[i]}') for i in range(len(u))]
        phyton_code += f"u = {u_completo}\n"
        phyton_code += f"f = {f}\n"
        phyton_code += f"fu = {fu}\n"
        phyton_code += f"hu = {hu}\n"
        phyton_code += f"fxw = {fxw}\n"
        phyton_code += f"hxw = {hxw}\n"
        phyton_code += "variables_locales = locals().copy()"
        with open(fullaffinename+'.py', 'w') as file:
            file.write(phyton_code)
        print(f"File {affine_modelname}.py has been created.")
       
        print('\n >>> Affine Conversion successful.\n')
        #Load affine model
        affine_model = importlib.import_module('models.affine_{}'.format(options.modelname))
        is_affine=1
#######################################################################################
    
    n  = len(affine_model.x)
    nu = len(affine_model.u)
    f1 = affine_model.fxw; #= 'f(x)' in eq. (2.1) in the paper
    f2 = affine_model.fu;  # = 'g(x)' in eq. (2.1) in the paper
    print(f"\nModel {options.modelname} has {n} states and {nu} inputs.\n")
########################################################################################
########################################################################################
       
# ======== Check for simbolic parameters and replace them ======================
    matrix_index = 1
    fcomp=[]
    while True:
        matrix_name = f'f{matrix_index}'
        if matrix_name in locals() or matrix_name in globals():
            matrix = locals().get(matrix_name, None)
            if matrix is not None:
                fcomp.append(Matrix(matrix))
        else:
            break
        matrix_index += 1
    vars = set()
    for element in fcomp:
        vars.update(element.free_symbols)
    import numpy as np
    k=np.array(list(vars))
    for i in range (0,len(model.x)):
        k=k[~np.isin(k,model.x)] 
    k=k.reshape(-1,1)
    random_values = {key[0]: np.random.randint(1, 11) for key in k}
    results_dict = {}
    for i, matrix in enumerate(fcomp, start=1):   
        result_matrix = matrix.subs(random_values)
        results_dict[f'f{i}'] = result_matrix
    f1 = results_dict['f1']
    f2= results_dict['f2']                
    print('Parameters replaced with random values')
    for key,value in random_values.items():
        print(f'{key} = {value}')
    t_end=time()
    elapsed_time=t_end-t_start
    print('(Preprocessing completed in {} seconds.)'.format(elapsed_time),end='\n\n')
    
###############Compute jacobians of f1 and f2##############################
    t_start=time()
    jac_f1=f1.jacobian(model.x)       
    jac_f2 = Matrix.zeros(n, n*nu)       
    for i in range(nu):
        jac_f2[:, i*n:(i+1)*n] = f2[:, i].jacobian(model.x)  
    t_end=time()
    elapsed_time=t_end-t_start
    print('(Compute jacobians completed in {} seconds.)'.format(elapsed_time))
#=========================================================================
#=== Compute the equilibrium point if one is not given ===================
    flattened_model_x = [item for sublist in model.x for item in sublist]
    if options.x0=='-':
        logic_eq=0
        sistema_homogeneo = [Eq(f, 0) for f in f1]
        soluciones = solve(sistema_homogeneo, flattened_model_x)    
        if isinstance(soluciones, dict):
            x0= [soluciones[var].evalf() for var in soluciones]
        else:    
            x0 = [valor for tupla in soluciones for valor in tupla]
    else:
        x0=options.x0
        if len(x0)!=len(model.x):
            print('Introduce a x0 vector with the proper dimensions')
            quit()
#=========================================================================
#=== Check the conditions indicated in 'opts': ===========================
    n = len(model.x)
    n_x0 = len(x0)/n
    if n_x0 > 1:        
        print('More than one equilibrium point')
        print('Condition check for the equilibrium point:')
        x0=x0[:n]
        print(f'    {x0}')
    else:
        if options.x0=='-':
            print('Condition check for the equilibrum point:')
            print(f'    {x0}')
        else:
            print('Condition check for the point:')
            print(f'    {x0}')
            print('If the point is not an equilibrium point, only the LARC results are valid.')
##################################################################
######## LC test #################################################
#######================================###########################
    if options.LCtest==1:
        tLC_start =time()  
        print("\n Computing the LINEARIZATION CONDITION (LC)... ") 
        A = f1.jacobian(model.x)
        B=f2
        variables_simbolicas = set()
        for elemento in A:
            variables_simbolicas.update(elemento.free_symbols)
        if variables_simbolicas:
            vars = A.free_symbols
            A_subs = A.subs(list(zip(vars, x0)))
        else:
            A_subs=A
        variables_simbolicas = set()
        for elemento in B:
            variables_simbolicas.update(elemento.free_symbols)
        if variables_simbolicas:
            vars = B.free_symbols
            B_subs = B.subs(list(zip(vars, x0)))
        else:
            B_subs=B           
        C_array=B_subs
        new_term=B_subs
        if options.numericLC==1:
            my_prime = nextprime(10**6)
            vars = A_subs.free_symbols.union(B.free_symbols)
            vars_num = [random.randint(0, my_prime) for _ in range(len(vars))]
            A_num = A_subs.subs(list(zip(vars, vars_num)))
            new_term_num = new_term.subs(list(zip(vars, vars_num)))
            C_array_num = C_array.subs(list(zip(vars, vars_num)))                   
           #--- Calculate the Controllability array numerically:           
            for i in range(0,n-1):               
                new_term_num=np.dot(A_num,new_term_num)               
                C_array_num = np.concatenate((C_array_num, new_term_num), axis=1)
           #--- Calculate rank:
            lin_rank = (Matrix(C_array_num)).rank()
        else: #--- Calculate the Controllability array symbolically:
            for i in range (0,n-1):
                new_term=np.dot(A_subs,new_term)
                C_array = np.concatenate((C_array, new_term), axis=1)
                lin_rank = (Matrix(C_array)).rank()
        print('\n rank(C(x)) =', lin_rank, end='')
        if lin_rank==n:
            print(f' => Model affine_{options.modelname} is short time locally controllable (STLC).')   
        else:
            print(' => The LC is not met');    
            print(f'=> it does not inform about the controllability of the model {options.modelname}\n');
        tLC_end =time()
        elapsed_time_LC=tLC_end-tLC_start
        print(f' (LC computation completed in {elapsed_time_LC} seconds)')

##################################################################
##################### ARC test####################################
###########===================================####################
# Loop that checks the Accesibility Rank Condition              
    if options.ARC==1:
        tARC_start=time()             
        print( '\n' +  'Computing the ACCESIBILITY RANK CONDITION (ARC)...\n Computing Lie bracket number')
        last_distr = np.zeros((n, sum(range(1, nu + 1))))
        def Lie_bracket(v1, v2, x, n, jac_v1):
            v1=Matrix(v1)
            v2=Matrix(v2)
            num_brcs = v1.shape[1] * v2.shape[1]
            brackets = Matrix()
            ibracket = 0
            for ivec2 in range(v2.shape[1]):
                for ivec1 in range(v1.shape[1]):
                    primer_factor_1=Matrix(v2[:, ivec2]).jacobian(model.x)
                    segundo_factor_1=Matrix(v1[:, ivec1])
                    primer_term=primer_factor_1*segundo_factor_1
                    primer_factor_2=Matrix(jac_v1[:, ivec1 * n:(ivec1 + 1) * n])
                    segundo_factor_2=Matrix(v2[:, ivec2]).reshape(n, 1)
                    segundo_term=primer_factor_2*segundo_factor_2
                    anadir=primer_term-segundo_term
                    brackets = brackets.row_join(Matrix(anadir))            
                    ibracket += 1                  
            return brackets
        # Assign values to first column up to 'nu' column of last_distr
        brackets = Lie_bracket(f1, f2, model.x, n, jac_f1)
        last_distr= np.hstack((last_distr, brackets))
        last_index=nu
        for j in range(1, nu):
            last_index = last_index + 1
            new_index = last_index + nu - j
            f2_j = f2[:, j]
            f2_rest = f2[:, j+1:]
            jac_f2_j = jac_f2[:, (j-1)*n:j*n]
            brackets = Lie_bracket(f2_j, f2_rest, x, n, jac_f2_j)
            last_distr[:, last_index:new_index] = brackets
        last_distr=np.sum(last_distr, axis=1)
        last_distr_expanded = np.expand_dims(last_distr, axis=1)
        Cx = np.concatenate((f1, f2, last_distr_expanded), axis=1)
        #Prov. test: evaluate at x=0 ***************************
        increaseLie=1
        iLie=1
        lasttime=0
        while increaseLie==1: #& lasttime< opts.maxtime:
            print(f' {iLie}... ')
            new_brackets = np.concatenate((Lie_bracket(f1, last_distr, model.x, n, jac_f1), Lie_bracket(f2, last_distr, model.x, n, jac_f2)),axis=1)                     
            #Construction of Brackets
            last_distr=np.transpose([(np.sum(new_brackets, axis=1))])
            Cx = np.concatenate((Cx, last_distr), axis=1)#ok<AGROW>
            if np.array_equal(last_distr, np.zeros((n, 1))):
                print(' => The ARC is not met')
                print(" => it does not inform about the accessibility of the model {}.\n".format(options.modelname))
                return
            
            if iLie >= math.ceil((n - 1) / nu):#If nu>1, it may not be necessary to calculate n-1 Lie brackets. The minimum we need to calculate is (n-1)/nu
                Cx_subs1=sym.Matrix(Cx)
                x0_dict = {var: val for var, val in zip(flattened_model_x, x0)}
                Cx_subs=Cx_subs1.subs(x0_dict)
                ctrl_rank=Cx_subs.rank()
                if ctrl_rank==n:
                    increaseLie=0;
            iLie= iLie+1
        if iLie-1 >= math.ceil((n - 1) / nu):
            if lasttime >= options.maxtime:
                print("\n    => More Lie brackets would be needed to see if the model is accessible.")
                print("\n    However, the maximum computation time allowed is reached.")
                print("\n    You can increase it by changing <<options.maxtime>> (currently options.maxtime = {})\n".format(options.maxtime))
                return
        print("\n rank(C(x)) = {}".format(ctrl_rank), end='')
        if ctrl_rank==n:
            print("   => Model {} is accessible.\n".format(options.modelname))
        else:
            if lasttime>= options.maxtime:
                print("\n    => More Lie brackets would be needed to see if the model is accessible.")
                print("\n    However, the maximum computation time allowed is reached.")
                print("\n    You can increase it by changing <<options.maxtime>> (currently options.maxtime = {})\n".format(options.maxtime))
            else:
                print("  => The ARC is not met")
                print(" => it does not inform about the accessibility of the model {}.\n".format(options.modelname))
        tARC_end=time()
        elapsed_time_ARC=tARC_end-tARC_start
        print('(ARC computation completed in {} seconds.)'.format(elapsed_time_ARC))
##################################################################
######## LARC test ###############################################
#######================================###########################
    if options.LARCtest==1 or options.GSC==1:
        print("\n Computing the LIE ALGEBRAIC RANK CONDITION (LARC)... ")  
        tLARC_start =time()
        def Lie_bracket(v1, v2, x, n, jac_v1):
            v1=Matrix(v1)
            v2=Matrix(v2)
            num_brcs = v1.shape[1] * v2.shape[1]
            brackets = Matrix()
            ibracket = 0
            for ivec2 in range(v2.shape[1]):
                for ivec1 in range(v1.shape[1]):
                    primer_factor_1=Matrix(v2[:, ivec2]).jacobian(model.x)
                    segundo_factor_1=Matrix(v1[:, ivec1])
                    primer_term=primer_factor_1*segundo_factor_1
                    primer_factor_2=Matrix(jac_v1[:, ivec1 * n:(ivec1 + 1) * n])
                    segundo_factor_2=Matrix(v2[:, ivec2]).reshape(n, 1)
                    segundo_term=primer_factor_2*segundo_factor_2
                    anadir=primer_term-segundo_term
                    brackets = brackets.row_join(Matrix(anadir))            
                    ibracket += 1                  
            return brackets
        #preprocessing
        Di=f1.row_join(f2)  #initial distribution
        x0_dict = {var: val for var, val in zip(flattened_model_x, x0)}
        Di_x0=Di.subs(x0_dict)
        Di_x0_simplified = Di_x0.applyfunc(simplify)
        Di_x0 = Di_x0_simplified.applyfunc(simplify)
        Di_rank=Di.rank()      #Rank of the initial distribution
        Di_x0_rank=Di_x0.rank()#Rank of the distribución in x0        
        Di_cols=Di.shape[1]    #Number of columns in the distribution
        Lie_index=[]                     #GSC index
        lie_index_added = False
        #lie_index_initialized = False
        last_brackets_index =eye(Di_cols)#Last brackets index of GSC
        last_brackets_index.col_del(0)
        last_brackets = f2               # First vectors to construct the brackets         
        last_brackets_numb = f2.shape[1]# Number of vectors to construct brackets
        # test to determine if the point is regular
        if Di_rank==Di_x0_rank:        #regularity of x0 in the initial distribution
            regular=1
        else:
            regular=0        
        #DISTRIBUTION COMPUTATION
        #Initialization
        increaseLie=1
        iLie=1
        lasttime=0
        print(f'Computing Lie bracket number {iLie}... ')
        base=Di
        base_index = list(range(1, base.shape[1] + 1))
        #Loop
        while increaseLie == 1 and regular==1 and lasttime < options.maxtime: #Iterations while not involutive and regular
            #===NEW BRACKETS COMPUTATION========#
            new_brackets = ([
                Lie_bracket(f1, last_brackets, flattened_model_x, n, jac_f1),
                Lie_bracket(f2, last_brackets, flattened_model_x, n, jac_f2)
           ])
            new_brackets=np.concatenate(new_brackets, axis=1)
            #indexes for the new brackets
            new_brackets_index=np.zeros((nu + 1, last_brackets_numb * (nu + 1)))
            #indexes for the new brackets
            new_brackets_index=np.zeros((nu + 1, last_brackets_numb * (nu + 1)))
            for i in range(nu + 1):
                aux_brackets_index = np.copy(last_brackets_index)
                filas, columnas = aux_brackets_index.shape
                eye_matrix = np.eye(1, last_brackets_numb)
                aux_brackets_index[i, :columnas] += eye_matrix.flatten()[:columnas]
                new_brackets_index[:, i * last_brackets_numb:(i + 1) * last_brackets_numb] = aux_brackets_index
            # number of new brackets
            new_brackets_numb= last_brackets_numb * (nu + 1)
            #remove null columns from the new brackets
            zero_cols_ind = []
            for icol in range(1,new_brackets.shape[1]+1): 
                if np.all(new_brackets[:, icol-1] == 0):  # Check for vectors of zeros
                    zero_cols_ind.append(icol)
            
            #Updates
            zero_cols_ind_adjusted = [index - 1 for index in zero_cols_ind]
            new_brackets = np.delete(new_brackets, zero_cols_ind_adjusted, axis=1)
            new_brackets_index = np.delete(new_brackets_index, zero_cols_ind_adjusted, axis=1)
            new_brackets_numb -= len(zero_cols_ind)
            #=======Compute new Distribution=====
            Di_np = np.array(Di)
            new_brackets_np = np.array(new_brackets)
            Dc_array = np.concatenate((Di_np, new_brackets_np), axis=1)
            Dc= Matrix(Dc_array)
            Dc=Dc.applyfunc(simplify)
            Dc=Dc.applyfunc(simplify)
            Dc=Dc.applyfunc(simplify)
            Dc_rank=(Dc).rank()
            Dc_x0 = np.array([[element.subs(x0_dict) for element in row] for row in Dc_array])
            Dc_x0_rank = (Matrix(Dc_x0)).rank()
            Dc_x0 = np.array([[element.subs(x0_dict) for element in row] for row in Dc_array])
            Dc_cols = Di_cols + new_brackets_numb
            #===VERIFICAITON OF DIFFERENT REQUIREMENTS====
            if new_brackets_numb == 0: # if no new non zero brackets STOP
                increaseLie=0
            #elif Dc_rank != Dc_x0_rank:     # if non regular point change loop
            #   regular=0
            else:
                if Dc_rank==n: #If full rank reached STOP
                    increaseLie=0
                elif Dc_rank==Di_rank: #If rank not increased STOP
                    increaseLie=0
                    #Updates
                    Dc=Di
                    Dc_x0=Di_x0
                    Dc_cols=Di_cols
                    new_brackets_numb=0
                    new_brackets_index=[]
                    new_brackets=[]
                if Dc_cols!=Dc_rank: #If colums dependent remove them
                    colum_depend_index=[]
                    Dc_check_old= Dc
                    Dc_check_old_rank=(Matrix(Dc).rank())
                    for icol in range(Di_cols + 1, Dc_cols):
                        Dc_check = Matrix(Dc)                                                        
                        Dc_check = np.delete(Dc_check,colum_depend_index,axis=1)
                        Dc_check_rank=(Matrix(Dc_check).rank())
                        if Dc_check_old_rank == Dc_check_rank:
                            colum_depend_index.extend([icol])
                            Dc_check_old = Dc.copy()
                            Dc_check_old = np.delete(Dc_check_old, colum_depend_index, axis=1)                         
                    #Updates
                    Dc_x0 = np.delete(Dc_x0, colum_depend_index, axis=1)

                    # Eliminar columnas especificadas de new_brackets y new_brackets_index
                    adjusted_columns = [idx - Di_cols for idx in colum_depend_index]
                    new_brackets = np.delete(new_brackets, adjusted_columns, axis=1)
                    new_brackets_index = np.delete(new_brackets_index, adjusted_columns, axis=1)
                    # Actualizar new_brackets_numb
                    new_brackets_numb -= len(colum_depend_index)
            #====UPDATE=====
            Di=Dc
            Di_x0=Dc_x0
            Di_rank=Dc_rank
            Di_cols=Di_cols+new_brackets_numb
            last_brackets=new_brackets
            if not lie_index_added:
                # Inicializar Lie_index con new_brackets_index si está vacío
                Lie_index = np.array(new_brackets_index)
                lie_index_added = True
            else:
                # Concatenar new_brackets_index como una columna adicional a Lie_index
                try:
                    new_brackets_index = (new_brackets_index[:, 0])
                    new_brackets_index = new_brackets_index.reshape(-1, 1)
                    Lie_index = np.concatenate((Lie_index, new_brackets_index), axis=1)
                except TypeError:
                    new_brackets_index = [item[0] for item in new_brackets_index if isinstance(item, tuple)]
                    n=4
                    tLARC_end =time()
                    elapsed_time_LARC=tLARC_end-tLARC_start
            last_brackets_index=new_brackets_index
            last_brackets_numb=new_brackets_numb
            iLie=iLie+1
            print(f' {iLie}... ', end='') 
        base=Di
        base_index=np.arange(1, base.shape[1] + 1)
        while increaseLie == 1 and regular==0 and lasttime < options.maxtime:
           # ======================NEW BRACKETS COMPUTATION===================== 
            new_brackets = ([
                Lie_bracket(f1, last_brackets, flattened_model_x, n, jac_f1),
                Lie_bracket(f2, last_brackets, flattened_model_x, n, jac_f2)
           ])
            new_brackets=np.concatenate(new_brackets, axis=1)
            #indexes for the new brackets
            new_brackets_index=np.zeros((nu + 1, last_brackets_numb * (nu + 1)))
            for i in range(nu + 1):
                aux_brackets_index = np.copy(last_brackets_index)
                aux_brackets_index[i, :] += np.ones(last_brackets_numb)
                new_brackets_index[:, i * last_brackets_numb:(i + 1) * last_brackets_numb] = aux_brackets_index
            # number of new brackets
            new_brackets_numb = np.array([[last_brackets_numb * (nu + 1)], [last_brackets_numb * (nu + 1)]])
            #======remove null colums from the new brackets=============
            zero_cols_ind = []
            for icol in range(new_brackets.shape[1]): 
                if np.all(new_brackets[:, icol] == 0):  # Check for vectors of zeros
                    zero_cols_ind.append(icol)
            #Updates
            new_brackets = np.delete(new_brackets, zero_cols_ind, axis=1)
            new_brackets_index = np.delete(new_brackets_index, zero_cols_ind, axis=1)
            new_brackets_numb = np.delete(new_brackets_numb, -len(zero_cols_ind), axis=0)
            zero_cols_ind = np.array(zero_cols_ind, dtype=np.int64)
            try:
                new_brackets_numb[:, 0] -= zero_cols_ind
            except ValueError:
                pass                    
            #------remove some of the dependent brackets----------
            #Initialization
            last_brackets = simplify(new_brackets)
            new_brackets = []
            last_brackets_numb = new_brackets_numb
            last_brackets_numb = int(np.ceil(last_brackets_numb[0, 0]))
            new_brackets_numb = 0
            last_brackets_index = new_brackets_index
            new_brackets_index = []
            Dc = Di
            #Find dependencies
            for i in range(1,last_brackets_numb+1):
                #Initializations
                logic=0
                base_cols= base.shape[1]
                base_index_dell=[]
                for j in range(1,base_cols+1):
                    #Increase gradually the base columns
                    aux=np.concatenate([base[:, -j:], last_brackets[:, i-1:i]], axis=1)
                    # Compute matrix of possible dependence
                    aux_matrix = Matrix(aux)
                    check = aux_matrix.nullspace()
                    if check:
                        check_end = check[0][:, -1]
                        try:
                        # Check for some smooth dependence
                            check_end_subs = check_end.subs(x, x0)
                            simplify(check_end_subs)
                            logic = 1  # Smooth dependence
                            break
                        except:
                            # Check for some smooth dependence to reduce base
                            index = np.where(check_end)[0]
                            if len(index) == 2:
                                auxx = aux.copy()
                                auxx[:, index[0]] = aux[:, index[1]]
                                auxx[:, index[1]] = aux[:, index[0]]
                                auxx_matrix = Matrix(auxx)
                                check_aux = auxx_matrix.nullspace()
                                try:
                                    # Reduction of the base
                                    check_aux_end_subs = check_aux[0][:, -1].subs(x, x0)
                                    simplify(check_aux_end_subs)
                                    base_index_dell.append(base_cols - j + index[0])
                                except:
                                    pass
                if logic==0: #dependence not found    
                    #updates
                    Dc= np.concatenate([Dc, last_brackets[:, i-1:i]], axis=1)  # Update distribution
                    base = np.delete(base, base_index_dell, axis=1)
                    base_index = [index for i, index in enumerate(base_index) if i not in base_index_dell]
                    base_index.append(Dc.shape[1])
                    base = np.concatenate([base, last_brackets[:, i-1:i]], axis=1)
                    new_brackets_numb += 1
                    new_brackets = new_brackets+last_brackets.tolist()
                    new_brackets_index =new_brackets_index+last_brackets_index.tolist()
            #Compute new distribution
            try:
                Dc = np.concatenate([Di, new_brackets], axis=1)             
            except ValueError:
                pass
            Dc_x0 = np.array([[element.subs(x0_dict) for element in row] for row in Dc])
            Dc_x0_rank = (Matrix(Dc_x0)).rank()
            Dc_cols = Di_cols + new_brackets_numb         
            #Verification of different requirements
            if new_brackets_numb == 0 or Dc_x0_rank == n: # If no new non zero brackets
                                                 # or full rank reached STOP 
                increaseLie = 0;
            #Update
            Di = Dc.copy()
            Di_x0 = Dc_x0.copy()
            Di_cols = Dc_cols
            last_brackets = new_brackets.copy()
            Lie_index += new_brackets_index
            last_brackets_index = new_brackets_index
            last_brackets_numb = new_brackets_numb
            iLie += 1
            print(f' {iLie}... ', end='')   
            tLARC_end =time()
            elapsed_time_LARC=tLARC_end-tLARC_start   

        
    #=======LARC display=======
    if options.LARCtest==1:
        print('\n dim(Delta_c(x)) = %d' % Dc_x0_rank, end='')
        if Dc_x0_rank==n: #reachable model
            tLARC_end =time()
            elapsed_time_LARC=tLARC_end-tLARC_start 
            print(f'   => The LARC is fulfilled => Model {options.modelname} is accessible.')
            print(f' (LARC computation completed in {elapsed_time_LARC} seconds)')
        else: #Unreachable model
            if elapsed_time_LARC> options.maxtime:
                print('\n    => More Lie brackets would be needed to see if the model is accessible.')
                print('\n    However, the maximum computation time allowed is reached.')
                print(f'\n    You can increase it by changing <<options.maxtime>> (currently options.maxtime = {options.maxtime})')
                print('\n    GSC can not be computed since LARC stopped before getting a conclusion')               
            else:  
                print(f'   => The LARC is not fulfilled => Model {options.modelname} is inaccessible.')
            print(f' (LARC computation completed in {elapsed_time_LARC} seconds)')
           
    #GSC display=====================================
    if options.GSC==1:
        tGSC_start=time()
        if Dc_x0_rank==n:
            Dc_x0_brackets=Dc_x0[:, nu+2:]
            try: 
                Lie_index = np.array(Lie_index)
            except ValueError:
                Lie_index=[[1,2],[1,1]]
                Lie_index = np.array(Lie_index)
            if not Lie_index.all():#GSC fulfilled
                print(f'   => The GSC is fulfilled => Model {options.modelname} is STLC.')
            else: # Bad brackets exist
                # Finding bad brackets
                bad_brackets_ind_aux = np.sum(np.abs(Dc_x0_brackets), axis=0)
                bad_brackets_ind_aux = np.nonzero(bad_brackets_ind_aux)[0]
                bad_brackets_ind = [[row[i] for i in bad_brackets_ind_aux] for row in Lie_index]
                bad_brackets_ind[0] = [np.mod(x, 2) for x in bad_brackets_ind[0]]
                for i in range(1, len(bad_brackets_ind)):
                    for j in range(len(bad_brackets_ind[i])):
                        bad_brackets_ind[i][j] = np.mod(1 + bad_brackets_ind[i][j], 2)
                bad_brackets_ind = np.sum(bad_brackets_ind, axis=0)
                bad_brackets_ind = nu + 1 - bad_brackets_ind
                bad_brackets_ind = bad_brackets_ind_aux[bad_brackets_ind == 0]
                if not bad_brackets_ind: #GSC fulfilled                                          
                    print()
                    print(f'   => The GSC is fulfilled => Model {options.modelname} is STLC.')
        tGSC_end =time()
        elapsed_time_GSC=tGSC_end-tGSC_start 
        print(f' (GSC check completed in {elapsed_time_GSC} seconds)')