function c2_write_aic_samples

dof_base = [24,32,32,32,40,40,40,48]+12+4;

genes = readtable('./data_input_all_samples/list_genes_all_samples.csv');




n = 80*3;

for n_gene = 1:1042




    temp = genes.Gene(n_gene);
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



    filename = ['./results_samples/' gene_name '.txt'];

    data = dlmread(filename);


    for line = 1:8

        chi = data(line,4);

        k = dof_base(line);


        aka = 2*k + chi + (2*k^2+2*k)/(n-k-1);

        dlmwrite(['./results_aic_samples/results_aic_samples_' gene_name  '.txt'],[line,chi,k,aka],'-append');

    end
end
end




