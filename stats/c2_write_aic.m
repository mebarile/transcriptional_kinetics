function c2_write_aic(pop1,pop2)
%% This code computes the akaike index per each gene

dof_base = [24,32,32,32,40,40,40,48]+8;

filePh = fopen('./data_input/list_genes.txt','r');
gene_names = textscan(filePh,'%s','delimiter','\n');
fclose(filePh);


for n = 1:1042

gene_name = gene_names{1}{n};


    filename = ['./results/gene_' gene_name '_' num2str(pop1) '_' num2str(pop2) '.txt'];

    data = dlmread(filename);


    for line = 1:8

        chi = data(line,4);

        k = dof_base(line);


        aka = 2*k + chi + (2*k^2+2*k)/(n-k-1);

        dlmwrite(['./results_aic/results_aic_' num2str(n1) '_' num2str(n2) '_' num2str(pop1) '_' num2str(pop2) '.txt'],[line,chi,k,aka],'-append');


    end
end


