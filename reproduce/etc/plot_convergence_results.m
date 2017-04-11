function plot_convergence_results(results)
fig = figure;
maxtime = 0;
colors = {[0.5 0.5 1],[0 0 0.7],[0 0 0.5],[0.7 0 0],[0.7 0 0],[0 0.7 0],[0.7 0 0],[0.7 0 0]};
markerstyles = {'+','^','o','s','*','+','+','+','^','o'};
for i = 1:length(results)
    thisresult = results{i};
    legnames{i} = thisresult.name;
    MSE = thisresult.stats.MSE;
    time = thisresult.stats.time.iter;
    cumtime = cumsum(time);
    semilogy(cumtime,MSE,'Color',colors{i},'LineStyle','-','LineWidth',1);
    hold all;
    maxtime = max(cumtime(end),maxtime);
end
hold off;
legend(legnames);
box on;
xlabel('CPU time (s)');
ylabel('NMSE');
axis tight
end
