function main_model3(n,reg,pop1,pop2,model)

addpath(genpath('./matlab_simulations/'))


rng(1)


filename = ['./data_output/parameters_' num2str(n) '_' num2str(model) '_model3_reg' num2str(reg) '_' num2str(pop1) '_' num2str(pop2) '.mat'];


load(['./data_input/data_gene_' num2str(n) '_model' num2str(model) '.mat'])



data_input = D.pop.mean;
data_input = data_input(:,[(pop1-1)*2+1,pop1*2,(pop2-1)*2+1,pop2*2]);

%% Definition of the Paramter Estimation Problem

parameters.name = {'aw1','aw2','aw3','aw4','aw5','aw6','aw7','aw8',...
    'bw1','bw2','bw3','bw4','bw5','bw6','bw7','bw8',...
    'gw1','gw2','gw3','gw4','gw5','gw6','gw7','gw8',...
    'bk1','bk2','bk3','bk4','bk5','bk6','bk7','bk8',...
    'Iwu','Iws','Iku','Iks','e1','e2','e3','e4'};

parameters.number = length(parameters.name);

modelfun = @ simulate_kinetics_model3;


parameters.min = [-11.873 * ones(32,1); log(0.1*(data_input(1,:)'+0.00001));log(0.01*(min(data_input)'+0.00001))];
parameters.max = [ones(32,1) * 6.7123;  log(3*(data_input(1,:)'+0.00001));  log(0.5*(max(data_input)'+0.00001))];

data1 = load(['data_output/parameters_' num2str(n) '_'  num2str(model) '_model1_reg' num2str(reg) '_' num2str(pop1) '_' num2str(pop2) '.mat']);
guess = data1.parameters.MS.par(:,1);
guess = guess([1:24,9:16,25:end]);

parameters.guess = guess;


% regularization hyperparameter
options.alpha = reg;


% Log-likelihood function
objectiveFunction = @(theta) llKin_model3(theta,modelfun,D,pop1,pop2,options);

%% Multi-start local optimization

% Options

optionsMultistart = PestoOptions();
optionsMultistart.obj_type = 'negative log-posterior';
optionsMultistart.n_starts = 3000;
optionsMultistart.comp_type = 'sequential';
optionsMultistart.mode = 'text';
optionsMultistart.proposal = 'uniform';
optionsMultistart.localOptimizerOptions.Display = 'iter';
optionsMultistart.localOptimizerOptions.Gradobj = 'on';
optionsMultistart.localOptimizerOptions.MaxIter = 6000;
optionsMultistart.localOptimizerOptions.MaxFunEvals = 12000;
% Optimization

parameters = getMultiStarts(parameters, objectiveFunction, optionsMultistart);


save(filename,'parameters','options','optionsMultistart')

end
