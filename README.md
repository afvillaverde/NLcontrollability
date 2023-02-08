# NLcontrollability
  MATLAB toolbox for analysing controllability and accessibility of nonlinear systems.


# To use NLcontrollability you need to:
  1. Download the code.
  2. Open a MATLAB session.
  3. Edit the 'run.m' file and run it.

  You can find an example in the 'run.m' file and more in 'ctrl_reach_unit_test_script.m'
    

# Options

  opts.LC        --> 1 if you want to check the Linearization Condition
  
                     0 if not
  opts.ARC       --> 1 if you want to check the Accessibility Rank Condition
                     0 if not
  opts.LARC      --> 1 if you want to check the Lie Algebraic Rank Condition
                     0 if not
  opts.GSC       --> 1 if you want to check Sussmann's General Sufficient Condition
                     0 if not
  'model_name'   --> Name of a .mat file placed in the 'models' folder
  opts.numericLC --> 0 if you want to check the Linearization Condition symbolically
                     1 for numeric computation
  opts.maxtime   --> max time for each test in seconds
  x0             --> Specific initial point. If no point is given it would try to compute equilibrium points.
      
               
# Entering the model
  NLcontrollability reads models stored as MATLAB MAT-files (.mat). Models are defined as follows:
  
  First, all parameters, states, and any other entities (such as inputs or known constants) appearing in the
  model must be defined as symbolic variables. Example:
    syms x1 x2 x3 x4 lambda u1
  Then we define the state variables, by creating a column vector named x. Example:
    x = [x1; x2; x3; x4];
  Similarly, we define the known input vector, u. If there are no inputs, enter a blank vector. Example:
    u = u1; % Or, for more than one input --> u = [u1; u2];
  The model ODEs (dx/dt=...) must also be entered as a column vector, called f, which must have the same
  length as the state vector x. Example:
    f = [u1;
		x1;
		x1^3;
		x3^2-lambda*x2^2*x1^4];
  Finally, save all the variables in a MAT-file. Example:
  save(‘MAPK’,‘x’,‘u’,‘f’);

  This model description format is the same used by the [STRIKE-GOLDD](https://github.com/afvillaverde/strike-goldd) toolbox. 
  Further details can be found in its [documentation](https://github.com/afvillaverde/strike-goldd/blob/master/STRIKE-GOLDD/doc/STRIKE-GOLDD_manual.pdf).

  
# Software contents
  **'run.m'**: file where the user enters the model and options. Running it executes the code.
  **'ctrl_analysis_MAIN.m':** main file; it calls the scripts with the accessibility and controllability tests.
  **'ctrl_LC.m'**: implementation of the Linearization Condition.
  **'ctrl_ARC.m'**: implementation of the Accessibility Rank Condition.
  **'ctrl_LARC_GSC.m'**: implementation of the Lie Algebraic Rank Condition and Sussmann's General Sufficient Condition.
  **'ctrl_reach_unit_tests_script.m'**: file with some examples ready to run.
  **'Lie_bracket.m'**: auxiliary function that computes the Lie bracket of two vectors w.r.t. a variable.
  **'frac_elem_sym.m'**: auxiliary function that changes a rational equation into its numerator and denominator.
  **'make_affine_known_u.m'**: function that converts a model into affine in inputs form


# Publications
  Publication of the methodology (pending)
