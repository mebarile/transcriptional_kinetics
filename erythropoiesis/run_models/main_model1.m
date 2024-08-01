function main_model1(n,reg,pop1,pop2)

addpath(genpath('./matlab_simulations/'))

rng(1)

filePh = fopen('./data_input/list_genes.txt','r');
gene_names = textscan(filePh,'%s','delimiter','\n');
fclose(filePh);


gene_name = gene_names{1}{n};

filename = ['./data_output/parameters_' gene_name '_model1_reg' num2str(reg) '_' num2str(pop1) '_' num2str(pop2) '.mat'];


load(['./data_input/data_gene_' gene_name '.mat'])

data_input = D.pop.mean;
data_input = data_input(:,[(pop1-1)*2+1,pop1*2,(pop2-1)*2+1,pop2*2]);

%% Definition of the Paramter Estimation Problem

parameters.name = {'aw1','aw2','aw3','aw4','aw5','aw6','aw7','aw8',...
    'bw1','bw2','bw3','bw4','bw5','bw6','bw7','bw8',...
    'gw1','gw2','gw3','gw4','gw5','gw6','gw7','gw8',...
    'Iwu','Iws','Iku','Iks','e1','e2','e3','e4'};

parameters.number = length(parameters.name);

modelfun = @ simulate_kinetics_model1;


parameters.min = [-11.873 * ones(24,1); log(0.1*(data_input(1,:)'+0.00001));log(0.01*(min(data_input)'+0.00001))];
parameters.max = [ones(24,1) * 6.7;  log(3*(data_input(1,:)'+0.00001));  log(0.1*(max(data_input)'+0.00001))];


% regularization hyperparameter
options.alpha = reg;


% Log-likelihood function
objectiveFunction = @(theta) llKin_model1(theta,modelfun,D,pop1,pop2,options);

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
% optionsMultistart.start_index = 6;


% Optimization

parameters = getMultiStarts(parameters, objectiveFunction, optionsMultistart);


save(filename,'parameters','options','optionsMultistart')
end
