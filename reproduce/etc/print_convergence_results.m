function print_convergence_results(results,tol)
%tol = 1e-4;
fprintf('Convergence criteria: MSE <= %2.2e\n',tol);
fprintf('algorithm  iter\t CPU time (s)\n',tol);
for i = 1:length(results)
    thisresult = results{i};
    MSE = thisresult.stats.MSE;
    time = thisresult.stats.time.iter;
    citer = 0;
    converge_flag = 0;
    for j=1:length(MSE);
        if( (MSE(j) < tol) && (citer == 0) )
            citer = j;
            converge_flag = 1;
        end
    end
    
    if converge_flag
        convergencetime = sum(time(1:citer));
    else
        convergencetime = Inf;
        citer = Inf;
    end
    fprintf('%9s: %d \t %.1f\n',thisresult.name,citer-1,convergencetime);
end
fprintf('\n');
end