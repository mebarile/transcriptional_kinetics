%% This code generates the input for modelling the erythropoietic dataset (Jack2 mutant), treating each replicate separately
% For each of the 1042 genes, the code creates a data structure D, which contains the following parts:
% D.pop.tw: the mean pseudotime of the cells in condition w (w is  the wild type)
% D.pop.tm: the mean pseudotime of the cells in condition m (m is the mutant)

% D.pop.t: the overall mean pseudotime, sorted
% D.pop.mean: mean spliced and unspliced counts in both conditions, for
% each sample
% D.pop.var: variance of the spliced and unspliced counts in both conditions




data_w_s1 = readtable('./data_input_all_samples/erythropoiesis_wt_spl_rep1.csv');
data_w_u1 = readtable('./data_input_all_samples_/erythropoiesis_wt_unspl_rep1.csv');

data_m_s1 = readtable('./data_input_all_samples/erythropoiesis_mut_spl_rep1.csv');
data_m_u1 = readtable('./data_input_all_samples/erythropoiesis_mut_unspl_rep1.csv');

data_w_s2 = readtable('./data_input_all_samples/erythropoiesis_wt_spl_rep2.csv');
data_w_u2 = readtable('./data_input_all_samples/erythropoiesis_wt_unspl_rep2.csv');

data_m_s2 = readtable('./data_input_all_samples/erythropoiesis_mut_spl_rep2.csv');
data_m_u2 = readtable('./data_input_all_samples/erythropoiesis_mut_unspl_rep2.csv');

data_w_s3 = readtable('./data_input_all_samples/erythropoiesis_wt_spl_rep3.csv');
data_w_u3 = readtable('./data_input_all_samples/erythropoiesis_wt_unspl_rep3.csv');

data_m_s3 = readtable('./data_input_all_samples/erythropoiesis_mut_spl_rep3.csv');
data_m_u3 = readtable('./data_input_all_samples/erythropoiesis_mut_unspl_rep3.csv');


genes = readtable('./data_input/list_genes.csv');
gene_names = genes.Gene;


%
data_w_s1 = sortrows(data_w_s1,'dpt_pseudotime');
data_m_s1 = sortrows(data_m_s1,'dpt_pseudotime');

data_w_u1 = sortrows(data_w_u1,'dpt_pseudotime');
data_m_u1 = sortrows(data_m_u1,'dpt_pseudotime');
%
data_w_s2 = sortrows(data_w_s2,'dpt_pseudotime');
data_m_s2 = sortrows(data_m_s2,'dpt_pseudotime');

data_w_u2 = sortrows(data_w_u2,'dpt_pseudotime');
data_m_u2 = sortrows(data_m_u2,'dpt_pseudotime');
%
data_w_s3 = sortrows(data_w_s3,'dpt_pseudotime');
data_m_s3 = sortrows(data_m_s3,'dpt_pseudotime');

data_w_u3 = sortrows(data_w_u3,'dpt_pseudotime');
data_m_u3 = sortrows(data_m_u3,'dpt_pseudotime');
%

scale_time = 1.42; % this just scales the dpt to fit in the same range as the model building.

nbins = 20;

for gene = 1:1042



    temp = genes.Gene(gene);
    gene_name = temp{1};

    % some gene names are changed when read from python to matlab


    if (length(gene_name) > 3) & ( (strcmp( gene_name(end-3:end-1),'Rik')) | (strcmp( gene_name(end-2:end),'Rik'))) &  ~(strcmp( gene_name(1),'C')) &  ~(strcmp( gene_name(1),'A')) &  ~(strcmp( gene_name(1),'D')) & ~(strcmp( gene_name(1),'B')) &  ~(strcmp( gene_name(1),'G'))  &  ~(strcmp( gene_name(1),'E')) &  ~(strcmp( gene_name(1),'F'))

        gene_name = ['x',gene_name]
    end

    if strcmp(gene_name,'Gt(ROSA)26Sor')

        gene_name = 'Gt_ROSA_26Sor'

    end


    gene_name(gene_name == '-') = '_';
    gene_name(gene_name == '.') = '_';


    out_name =  ['./data_input_all_samples/data_gene_' gene_name '.mat'];


    measured_w1 = [data_w_u1{:,gene_name},data_w_s1{:,gene_name},data_w_s1{:,'dpt_pseudotime'}];
    measured_m1 = [data_m_u1{:,gene_name},data_m_s1{:,gene_name},data_m_s1{:,'dpt_pseudotime'}];

    measured_w2 = [data_w_u2{:,gene_name},data_w_s2{:,gene_name},data_w_s2{:,'dpt_pseudotime'}];
    measured_m2 = [data_m_u2{:,gene_name},data_m_s2{:,gene_name},data_m_s2{:,'dpt_pseudotime'}];

    measured_w3 = [data_w_u3{:,gene_name},data_w_s3{:,gene_name},data_w_s3{:,'dpt_pseudotime'}];
    measured_m3 = [data_m_u3{:,gene_name},data_m_s3{:,gene_name},data_m_s3{:,'dpt_pseudotime'}];

    %% group

    remove = mod(size(measured_w1,1),nbins)+20;
    measured_w1(1:remove,:) = [];

    n_w1 = size(measured_w1,1);

    %

    remove = mod(size(measured_m1,1),nbins)+40;
    measured_m1(1:remove,:) = [];

    n_m1 = size(measured_m1,1);


    %% compute means w1

    spliced_w1 = reshape(measured_w1(:,2),n_w1/nbins,nbins);
    unspliced_w1 = reshape(measured_w1(:,1),n_w1/nbins,nbins);

    mean_spliced_w1 = mean(spliced_w1);
    mean_unspliced_w1 = mean(unspliced_w1);

    var_spliced_w1 = var(spliced_w1);
    var_unspliced_w1 = var(unspliced_w1);

    time_w1 = reshape(measured_w1(:,end),n_w1/nbins,nbins);
    mean_time_w1 = mean(time_w1);


    mean_time_w1 = mean_time_w1 * ww;



    %% compute means m1

    spliced_m1 = reshape(measured_m1(:,2),n_m1/nbins,nbins);
    unspliced_m1 = reshape(measured_m1(:,1),n_m1/nbins,nbins);

    mean_spliced_m1 = mean(spliced_m1);
    mean_unspliced_m1 = mean(unspliced_m1);

    var_spliced_m1 = var(spliced_m1);
    var_unspliced_m1 = var(unspliced_m1);

    time_m1 = reshape(measured_m1(:,end),n_m1/nbins,nbins);
    mean_time_m1 = mean(time_m1);

    mean_time_m1 = mean_time_m1 * ww;



    %% group

    remove = mod(size(measured_w2,1),nbins)+20;
    measured_w2(1:remove,:) = [];

    n_w2 = size(measured_w2,1);

    %

    remove = mod(size(measured_m2,1),nbins)+40;
    measured_m2(1:remove,:) = [];

    n_m2 = size(measured_m2,1);


    %% compute means w2

    spliced_w2 = reshape(measured_w2(:,2),n_w2/nbins,nbins);
    unspliced_w2 = reshape(measured_w2(:,1),n_w2/nbins,nbins);

    mean_spliced_w2 = mean(spliced_w2);
    mean_unspliced_w2 = mean(unspliced_w2);

    var_spliced_w2 = var(spliced_w2);
    var_unspliced_w2 = var(unspliced_w2);

    time_w2 = reshape(measured_w2(:,end),n_w2/nbins,nbins);
    mean_time_w2 = mean(time_w2);

    mean_time_w2 = mean_time_w2 * ww;

    %% compute means m2

    spliced_m2 = reshape(measured_m2(:,2),n_m2/nbins,nbins);
    unspliced_m2 = reshape(measured_m2(:,1),n_m2/nbins,nbins);

    mean_spliced_m2 = mean(spliced_m2);
    mean_unspliced_m2 = mean(unspliced_m2);

    var_spliced_m2 = var(spliced_m2);
    var_unspliced_m2 = var(unspliced_m2);

    time_m2 = reshape(measured_m2(:,end),n_m2/nbins,nbins);
    mean_time_m2 = mean(time_m2);

    mean_time_m2 = mean_time_m2 * ww;


    %% group

    remove = mod(size(measured_w3,1),nbins)+20;
    measured_w3(1:remove,:) = [];

    n_w3 = size(measured_w3,1);

    %

    remove = mod(size(measured_m3,1),nbins)+40;
    measured_m3(1:remove,:) = [];

    n_m3 = size(measured_m3,1);


    %% compute means w3

    spliced_w3 = reshape(measured_w3(:,2),n_w3/nbins,nbins);
    unspliced_w3 = reshape(measured_w3(:,1),n_w3/nbins,nbins);

    mean_spliced_w3 = mean(spliced_w3);
    mean_unspliced_w3 = mean(unspliced_w3);

    var_spliced_w3 = var(spliced_w3);
    var_unspliced_w3 = var(unspliced_w3);

    time_w3 = reshape(measured_w3(:,end),n_w3/nbins,nbins);
    mean_time_w3 = mean(time_w3);

    mean_time_w3 = mean_time_w3 * ww;

    %% compute means m3

    spliced_m3 = reshape(measured_m3(:,2),n_m3/nbins,nbins);
    unspliced_m3 = reshape(measured_m3(:,1),n_m3/nbins,nbins);

    mean_spliced_m3 = mean(spliced_m3);
    mean_unspliced_m3 = mean(unspliced_m3);

    var_spliced_m3 = var(spliced_m3);
    var_unspliced_m3 = var(unspliced_m3);

    time_m3 = reshape(measured_m3(:,end),n_m3/nbins,nbins);
    mean_time_m3 = mean(time_m3);

    mean_time_m3 = mean_time_m3 * ww;

    %% for fit

    measured_w = [mean_unspliced_w1', mean_spliced_w1', mean_unspliced_w2', mean_spliced_w2', mean_unspliced_w3', mean_spliced_w3'];
    var_w      = [var_unspliced_w1' , var_spliced_w1' , var_unspliced_w2' , var_spliced_w2' , var_unspliced_w3' , var_spliced_w3' ];

    measured_m = [mean_unspliced_m1', mean_spliced_m1', mean_unspliced_m2', mean_spliced_m2', mean_unspliced_m3', mean_spliced_m3'];
    var_m      = [var_unspliced_m1' , var_spliced_m1' , var_spliced_m2'   , var_spliced_m2' , var_spliced_m3'   , var_spliced_m3' ];




    %%


    D.pop.tw = [mean_time_w1' - mean_time_w1(1), mean_time_w2' - mean_time_w2(1), mean_time_w3' - mean_time_w3(1)];
    D.pop.tm = [mean_time_m1' - mean_time_m1(1), mean_time_m2' - mean_time_m2(1), mean_time_m3' - mean_time_m3(1)];

    D.pop.t = sort([mean_time_w1-mean_time_w1(1), mean_time_w2-mean_time_w2(1), mean_time_w3-mean_time_w3(1), ...
        mean_time_m1-mean_time_m1(1), mean_time_m2-mean_time_m2(1), mean_time_m3-mean_time_m3(1)]);

    D.pop.mean = [measured_w,measured_m];
    D.pop.var = [var_w,var_m];


    %'saving'
    save(out_name,'D');
end
