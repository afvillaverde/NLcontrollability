# NLcontrollability_py
  Python scripts for analysing controllability and accessibility of nonlinear systems.


# To use NLcontrollability you need to:
  1. Download the code.
  2. Put all the files in the same folder.
  3. Open an application capable of executing Python code.
  4. Edit the 'options.py' file and run it.
  5. Run the file 'runPY'.

# Options
 -**modelname**      --> Write the name of a model located in the 'models' folder.
   - **LCtest**        --> 1 if you want to check the Linearization Condition, 0 if not
   - **numericLC** -->0 if you want to check the Linearization Condition symbolically, 1 for numeric computation
   - **ARC**       --> 1 if you want to check the Accessibility Rank Condition, 0 if not
   - **LARCtest**      --> 1 if you want to check the Lie Algebraic Rank Condition, 0 if not
   - **GSC**       --> 1 if you want to check Sussmann's General Sufficient Condition, 0 if not
   
   - **x0**             --> Specific initial point. If no point is given (x0=¨-¨ it would try to compute equilibrium points.
      
               
# Entering the model
  NLcontrollability reads models stored as PYTHON Py-files (.py). Models are defined as follows:
  
  First, all parameters, states, and any other entities (such as inputs or known constants) appearing in the
  model must be defined as symbolic variables. Example:
  
	x1, x2, x3, u1, p1 = sym.symbols('x1 x2 x3 u1 p1')
	
  Then we define the state variables, by creating a column vector named x. Example:
  
	x = [[x1], [x2], [x3]]

  Similarly, we define the known input vector, u. If there are no inputs, enter a blank vector. Example:
  
	u = [u1]; % Or, for more than one input --> u = [u1, u2];
	
Define the differential equations (ODEs) for the model as a column vector f. This vector should have the same length as x and represent dx/dt equations. For example:
  
f = [[u1],
     [x1],
     [x1***3],
     [x3 - p1*x2]]

		
  Finally, store all variables locally for debugging or further use.
  
  This model description format is the same used by the [STRIKE-GOLDD](https://github.com/afvillaverde/strike-goldd) toolbox. 
  Further details can be found in its [documentation](https://github.com/afvillaverde/strike-goldd/blob/master/STRIKE-GOLDD/doc/STRIKE-GOLDD_manual.pdf).

  
# Software contents
  - **`run.py`**: file where the user enters the model. Running it executes the code.
  - **`options.py`**: file where the user selects the analysis to conduct.
  - **`nl_controllability.py`:** main file.;


