function c1_analyse_genes_samples
%% this code puts together the results for each gene

addpath(genpath('./matlab_simulations'))
genes = readtable('./data_input_all_samples/list_genes_all_samples.csv');


for n = 1:1042


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

    name = gene_name;




    for model = 1:8

        filename = ['./data_output_all_samples/parameters_' name '_model' num2str(model) '_all_samples_reg10.mat'];
        data = load(['./data_input_all_samples/data_gene_' name '.mat']);


        data_input = data.D.pop.mean;


        t = data.D.pop.t;

        err_regulation = 0.1 * max(data_input);
        err_regulation  = repmat(err_regulation,20,1);

        para = load(filename);




        name_m =  ['simulate_kinetics_model' num2str(model) '_all_samples'];

        modelfun = str2func(name_m);


        theta = para.parameters.MS.par(:,1);


        npara = length(theta)-4;

        [~,~,~,N,~,us] = modelfun(t,theta,[],[],[]);



        ind_w1(1,1) = 1;
        ind_k1(1,1) = 1;

        ind_w2(1,1) = 1;
        ind_k2(1,1) = 1;

        ind_w3(1,1) = 1;
        ind_k3(1,1) = 1;



        time_all = [data.D.pop.tw,data.D.pop.tm];



        for iHelp = 2:20

            ind_w1(iHelp,1) = find(data.D.pop.t == time_all(iHelp,1));
            ind_k1(iHelp,1) = find(data.D.pop.t == time_all(iHelp,4));

            ind_w2(iHelp,1) = find(data.D.pop.t == time_all(iHelp,2));
            ind_k2(iHelp,1) = find(data.D.pop.t == time_all(iHelp,5));

            ind_w3(iHelp,1) = find(data.D.pop.t == time_all(iHelp,3));
            ind_k3(iHelp,1) = find(data.D.pop.t == time_all(iHelp,6));

        end


        N = [N(ind_w1,1:2),N(ind_w2,3:4),N(ind_w3,5:6),N(ind_k1,7:8),N(ind_k2,9:10),N(ind_k3,11:12)];


        negLogL = 0;

        diff =0;

        for it = 1:20
            for col = 1:12

                err_idx1 = floor((col-1)/6);
                err_idx2 = 2-mod(col,2);

                err_idx = 2*err_idx1+err_idx2;


                erri = exp(theta(npara+err_idx));


                negLogL = negLogL + 0.5*log(2*pi*erri^2) + ...
                    0.5*(data_input(it,col)-N(it,col)).^2/...
                    (erri^2);


                diff = diff + (data_input(it,col)-N(it,col)).^2/(err_regulation(it,col)^2);

            end
        end

        dlmwrite(['./results_samples/' name  '.txt'],[model,para.parameters.MS.logPost(1),negLogL,diff],'-append')

    end
end
end


