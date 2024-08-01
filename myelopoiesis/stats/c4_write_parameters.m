function c4_write_parameters(pop1,pop2)
%% This code writes the parameters estimated for each gene

filePh = fopen('./data_input/list_genes.txt','r');
gene_names = textscan(filePh,'%s','delimiter','\n');
fclose(filePh);

for n = 1:1228
    gene_name = gene_names{1}{n};



    for model = 1:8

        filename = ['./data_output/parameters_' gene_name '_model' num2str(model) '_reg10_' num2str(pop1) '_' num2str(pop2) '.mat'];
        data = load(filename);
        para = data.parameters.MS.par(:,1)';


        dlmwrite(['./results/fit_parameters.txt',[n,model,para],'-append')

    end

end

end

