%% This code generates the input for modelling the erythropoietic dataset (Jack2 mutant)
% For each of the 1042 genes, the code creates a data structure D, which contains the following parts:
% D.pop.tw: the mean pseudotime of the cells in condition w (w is  the wild type)
% D.pop.tm: the mean pseudotime of the cells in condition m (m is the mutant)

% D.pop.t: the overall mean pseudotime, sorted
% D.pop.mean: mean spliced and unspliced counts in both conditions
% D.pop.var: variance of the spliced and unspliced counts in both conditions


data_w_s = readtable('./data_input/jak2_wt_ery_spl.csv');
data_w_u = readtable('./data_input/jak2_wt_ery_unspl.csv');

data_m_s = readtable('./data_input/jak2_mut_ery_spl.csv');
data_m_u = readtable('./data_input/jak2_mut_ery_unspl.csv');

genes = readtable('./data_input/list_genes.csv');
gene_names = genes.Gene;


data_w_s = sortrows(data_w_s,'dpt_pseudotime');
data_m_s = sortrows(data_m_s,'dpt_pseudotime');

data_w_u = sortrows(data_w_u,'dpt_pseudotime');
data_m_u = sortrows(data_m_u,'dpt_pseudotime');

dlmwrite('./data_output/dpt_w.txt',data_w_s{:,'pseudotime'})
dlmwrite('./data_output/dpt_m.txt',data_m_s{:,'pseudotime'})


%%

scale_time = 1.42; % this just scales the dpt to fit in the same range as the model building. 

nbins = 20;


for gene = 1:1042



    temp = genes.Gene(gene);
    gene_name = temp{1};

    % some gene names are changed when read from python to matlab

    if (length(gene_name) > 3) & ( (strcmp( gene_name(end-3:end-1),'Rik')) | (strcmp( gene_name(end-2:end),'Rik'))) &  ~(strcmp( gene_name(1),'C')) &  ~(strcmp( gene_name(1),'A')) &  ~(strcmp( gene_name(1),'B')) &  ~(strcmp( gene_name(1),'G'))  &  ~(strcmp( gene_name(1),'F'))

        gene_name = ['x',gene_name];
    end

    gene_name(gene_name == '-') = '_';
    gene_name(gene_name == '.') = '_';


    out_name =  ['./data_input_2033/data_gene_' gene_name '.mat'];

    measured_w = [data_w_u{:,gene_name},data_w_s{:,gene_name},data_w_s{:,'dpt_pseudotime'}];
    measured_m = [data_m_u{:,gene_name},data_m_s{:,gene_name},data_m_s{:,'dpt_pseudotime'}];


    %% group

    remove = mod(size(measured_w,1),nbins)+20;
    measured_w(1:remove,:) = [];

    n_w = size(measured_w,1);

   

    remove = mod(size(measured_m,1),nbins)+40;
    measured_m(1:remove,:) = [];

    n_m = size(measured_m,1);


    %% compute means w

    spliced_w = reshape(measured_w(:,2),n_w/nbins,nbins);
    unspliced_w = reshape(measured_w(:,1),n_w/nbins,nbins);

    mean_spliced_w = mean(spliced_w);
    mean_unspliced_w = mean(unspliced_w);

    var_spliced_w = var(spliced_w);
    var_unspliced_w = var(unspliced_w);

    time_w = reshape(measured_w(:,end),n_w/nbins,nbins);
    mean_time_w = mean(time_w);


    mean_time_w = mean_time_w * scale_time;



    %% compute means m

    spliced_m = reshape(measured_m(:,2),n_m/nbins,nbins);
    unspliced_m = reshape(measured_m(:,1),n_m/nbins,nbins);

    mean_spliced_m = mean(spliced_m);
    mean_unspliced_m = mean(unspliced_m);

    var_spliced_m = var(spliced_m);
    var_unspliced_m = var(unspliced_m);

    time_m = reshape(measured_m(:,end),n_m/nbins,nbins);
    mean_time_m = mean(time_m);

    mean_time_m = mean_time_m * scale_time;

    %% for fit

    measured_w = [mean_unspliced_w',mean_spliced_w'];
    var_w = [var_unspliced_w',var_spliced_w'];



    measured_m = [mean_unspliced_m',mean_spliced_m'];
    var_m = [var_unspliced_m',var_spliced_m'];


    D.pop.tw = mean_time_w'-mean_time_w(1);
    D.pop.tm = mean_time_m'-mean_time_m(1);



    D.pop.t = sort([mean_time_w-mean_time_w(1),mean_time_m-mean_time_m(1)]);

    D.pop.mean = [measured_w,measured_m];
    D.pop.var = [var_w,var_m];
    
    %%'saving'
    save(out_name,'D');
end
