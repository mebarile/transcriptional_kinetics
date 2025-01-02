function c3_write_statistics_aic_samples


genes = readtable('./data_input_all_samples/list_genes_all_samples.csv');


for n = 1:1042


    temp = genes.Gene(n);
    gene_name = temp{1};

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



    filename = ['./results_aic_samples/res_aka_samples_' gene_name '.txt'];

    data = dlmread(filename);



    chi = data(1:8,end);

    [~, min_chi] = min(chi);


    dlmwrite(['./results_samples/statistic_aic_samples.txt'],[n,min_chi],'-append');

end

end

