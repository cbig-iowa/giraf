function results = run_experiments(x0,b,A,At,sampmask,exp)
results = cell(size(exp));
num_exp = length(results);
for i=1:num_exp
opt = exp{i}.opt_handle;
[x,stats] = opt(x0,b,A,At,sampmask,exp{i}.param,exp{i}.settings);
results{i}.name = exp{i}.name;
results{i}.x = x;
results{i}.stats = stats;
end
end