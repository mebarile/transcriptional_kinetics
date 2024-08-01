function c3_write_statistics_aic(pop1,pop2)
%% This code finds the best model for each gene



filePh = fopen('./data_input/list_genes.txt','r');
gene_names = textscan(filePh,'%s','delimiter','\n');
fclose(filePh);


for n = 1:1042
    gene_name = gene_names{1}{n};


    filename = ['./results/results_aic_' gene_name '_' num2str(pop1) '_' num2str(pop2) '.txt'];


    data = dlmread(filename);



    chi = data(1:8,end);

    [~, min_chi] = min(chi);

    if sum(isinf(chi))


        dlmwrite(['statistics_aic_' num2str(pop1) '_' num2str(pop2) '.txt'],[n,0],'-append');

    else
        dlmwrite(['./results/statistics_aic_' num2str(pop1) '_' num2str(pop2) '.txt'],[n,min_chi],'-append');
    end

end

end


