function main_model8_samples(n,reg)
addpath(genpath('./matlab_simulations/'))

rng(1)


genes = readtable('./data_input_all_samples/list_genes_all_samples.csv');



temp = genes.Gene(n);
gene_name = temp{1}

if (length(gene_name) > 3) & ( (strcmp( gene_name(end-3:end-1),'Rik')) | (strcmp( gene_name(end-2:end),'Rik'))) &  ~(strcmp( gene_name(1),'C')) & ...
        ~(strcmp( gene_name(1),'A')) &  ~(strcmp( gene_name(1),'D')) & ~(strcmp( gene_name(1),'B')) ...
        &  ~(strcmp( gene_name(1),'G'))  &  ~(strcmp( gene_name(1),'E')) &  ~(strcmp( gene_name(1),'F'))

    gene_name = ['x',gene_name]
end

if strcmp(gene_name,'Gt(ROSA)26Sor')

    gene_name = 'Gt_ROSA_26Sor'

end

gene_name(gene_name == '-') = '_';
gene_name(gene_name == '.') = '_';


filename = ['./data_output_all_samples/parameters_' gene_name '_model8_all_samples_reg' num2str(reg) '.mat'];


load(['./data_input_all_samples/data_gene_' gene_name '.mat'])

data_input = D.pop.mean;

%% Definition of the Paramter Estimation Problem

parameters.name = {'aw1','aw2','aw3','aw4','aw5','aw6','aw7','aw8',...
    'bw1','bw2','bw3','bw4','bw5','bw6','bw7','bw8',...
    'gw1','gw2','gw3','gw4','gw5','gw6','gw7','gw8',...
    'ak1','ak2','ak3','ak4','ak5','ak6','ak7','ak8',...
    'bk1','bk2','bk3','bk4','bk5','bk6','bk7','bk8',...
    'gk1','gk2','gk3','gk4','gk5','gk6','gk7','gk8',...
    'Iwu1','Iws1','Iwu2','Iws2','Iwu3','Iws3','Iku1','Iks1','Iku2','Iks2','Iku3','Iks3','e1','e2','e3','e4'};

parameters.number = length(parameters.name);

modelfun = @ simulate_kinetics_model8_all_samples;

parameters.min = [-11.873 * ones(48,1); log(0.1*(data_input(1,:)'+0.00001));log(0.01*(min(data_input(:,[1,3,5]),[],'all')'+0.00001));...
    log(0.01*(min(data_input(:,[2,4,6]),[],'all')'+0.00001));log(0.01*(min(data_input(:,[7,9,11]),[],'all')'+0.00001));log(0.01*(min(data_input(:,[8,10,12]),[],'all')'+0.00001))];

parameters.max = [ones(48,1) * 6.7123;  log(3*(data_input(1,:)'+0.00001));  log(0.1*(max(data_input(:,[1,3,5]),[],'all')'+0.00001));...
    log(0.1*(max(data_input(:,[2,4,6]),[],'all')'+0.00001));  log(0.1*(max(data_input(:,[7,9,11]),[],'all')'+0.00001)) ;log(0.1*(max(data_input(:,[8,10,12]),[],'all')'+0.00001))];


data5 = load(['./data_output_all_samples/parameters_' gene_name '_model5_all_samples_reg' num2str(reg) '.mat']);
b5 = data5.parameters.MS.logPost(1);


data6 = load(['./data_output_all_samples/parameters_' gene_name '_model6_all_samples_reg' num2str(reg) '.mat']);
b6 = data6.parameters.MS.logPost(1);


data7 = load(['./data_output_all_samples/parameters_' gene_name '_model7_all_samples_reg' num2str(reg) '.mat']);
b7 = data7.parameters.MS.logPost(1);

if (b5>=b6) & (b5>=b7)
    guess = data5.parameters.MS.par(:,1);
    guess = guess([1:40,17:24,41:end]);

elseif (b6>=b7) & (b6>=b5)
    guess = data6.parameters.MS.par(:,1);
    guess = guess([1:32,9:16,33:end]);

else

    guess = data7.parameters.MS.par(:,1);
    guess = guess([1:24,1:8,25:end]);
end

parameters.guess = guess;

% regularization hyperparameter
options.alpha = reg;


% Log-likelihood function
objectiveFunction = @(theta) llKin_model8_all_samples(theta,modelfun,D,options);

%% Multi-start local optimization

% Options

optionsMultistart = PestoOptions();
optionsMultistart.obj_type = 'negative log-posterior';
optionsMultistart.n_starts = 200;
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
end