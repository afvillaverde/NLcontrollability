%=========================================================================%
%= Function that checks the Accesibility Rank Condition ===============%
function ctrl_ARC(modelname,opts,n,nu,f1,f2,x,x0)

tARC = tic;
fprintf('\n Computing the ACCESIBILITY RANK CONDITION (ARC)...\n Computing Lie bracket number');
last_distr  =  sum(Lie_bracket(f1,f2,x,n),2) ;
Cx(:,1:nu+2)  = [f1,f2,last_distr];
%****** Prov. test: evaluate at x = 0 *************************************
increaseLie = 1;
iLie        = 1;
lasttime = 0;
while increaseLie == 1 && lasttime < opts.maxtime
    fprintf(' %d... ',iLie)
    new_brackets = [Lie_bracket(f1,last_distr,x,n), ... % Construction 
                    Lie_bracket(f2,last_distr,x,n)];    % of brackets
    last_distr = sum(new_brackets,2);
    Cx= [Cx,last_distr]; %#ok<AGROW> 
    if last_distr == zeros(n,1)
        fprintf('   => The ARC is not met');    
        fprintf(' => it does not inform about the accessibility of the model %s.\n', modelname);
        return
    end
    if iLie>=ceil((n-1)/nu) % If nu>1, it may not be necessary to calculate n-1 Lie brackets. The minimum we need to calculate is (n-1)/nu
        ctrl_rank   = rank(subs(Cx,x,x0));
        if ctrl_rank == n 
            increaseLie = 0;
        end
    end
    iLie = iLie+1    ;   
    lasttime = toc(tARC);
end      

fprintf('\n rank(C(x)) = %d',ctrl_rank);
if ctrl_rank == n
    fprintf('   => Model %s is accessible.\n', modelname);    
else
    if lasttime >= opts.maxtime
        fprintf('\n    => More Lie brackets would be needed to see if the model is accessible.');
        fprintf('\n    However, the maximum computation time allowed is reached.');
        fprintf('\n    You can increase it by changing <<opts.maxtime>> (currently opts.maxtime = %d)\n',opts.maxtime);
    else
        fprintf('   => The ARC is not met');    
        fprintf(' => it does not inform about the accessibility of the model %s.\n', modelname);   
    end
    fprintf(' (ARC computation completed in %d seconds.)\n ',toc(tARC));
end
end