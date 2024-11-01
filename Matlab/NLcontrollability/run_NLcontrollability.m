% running example for NLcontrolability

clear; clc;
addpath('models')

% ======================== Explanation of variables =======================
% 
% ctrl_analysis_MAIN([opts.LC,opts.ARC,opts.LARC,opts.GSC],'model_name',...
%                     opts.numericLC, opts.maxtime, x0(optional));
%     opts.LC        --> 1 if you want to check the Linearization Condition
%                        0 if not
%     opts.ARC       --> 1 if you want to check the Accesibility Rank 
%                          Condition
%                        0 if not
%     opts.LARC      --> 1 if you want to check the Lie Algebraic Rank 
%                          Condition
%                        0 if not
%     opts.GSC       --> 1 if you want to check Sussmanns General 
%                          Sufficient Condition
%                        0 if not
%     'model_name'   --> Name of a .mat file placed in the 'models' folder
%     opts.numericLC --> 0 if you want to check the Linearization
%                          Condition symbolicaly
%                        1 for numeric computation
%     opts.maxtime   --> max time for the each test in seconds
%     x0             --> Specific initial point. If no  point is given it 
%                        would try to compute equilibrium points.

% ================================ Example ================================
%
% You can rewrite the following line to run the analisys on your model,
% more example are given in ctrl_reach_unit_tests_script.m
ctrl_analysis_MAIN([0,1,1,0],'Fujita_try_xp', 0, 600,ones(9,1)); 

