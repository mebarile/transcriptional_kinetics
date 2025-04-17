function plot_all_models_best
%% read data

genes = readtable('./data_input/list_genes.csv');


pop1 = 1;
pop2 = 2;

stat = dlmread(['./results/statistics_aic_' num2str(pop1) '_' num2str(pop2) '.txt']);


for n = 490

    temp = genes.Gene(n);

    gene_name = temp{1}

    if (length(gene_name) > 3) & ( (strcmp( gene_name(end-3:end-1),'Rik')) | (strcmp( gene_name(end-2:end),'Rik'))) &  ~(strcmp( gene_name(1),'C')) & ...
            ~(strcmp( gene_name(1),'A')) &  ~(strcmp( gene_name(1),'B')) &  ~(strcmp( gene_name(1),'E'))  & ~(strcmp( gene_name(1),'G'))  &  ~(strcmp( gene_name(1),'F'))

        gene_name = ['x',gene_name]
    end
    best_model = stat(n,2);


    data = load(['./data_input/data_gene_' gene_name '.mat']);

    all_para  = dlmread(['./results/parameters_aka_' num2str(pop1) '_' num2str(pop2) '.txt']);


    time_all = [data.D.pop.tw,data.D.pop.tm];


    mean_time_wt = time_all(:,pop1);
    mean_time_ko = time_all(:,pop2);

    %% start fit



    data_input = data.D.pop.mean;
    measured = data_input(:,[(pop1-1)*2+1,pop1*2,(pop2-1)*2+1,pop2*2]);


    % plot_time = vec_time(1):0.001:vec_time(end);
    in = [29,37,37,37,45,45,45,53];

    %%
    figure(11)
    clf


    for model = best_model


        theta = all_para((n-1)*8+model,3:end);

        err_t = exp(theta(in(model):in(model)+3));
        error = repmat(err_t,20,1);



        name_m =  ['simulate_kinetics_4cond_model' num2str(model) ];

        modelfun = str2func(name_m);



        %% model


        [~,~,~,N] = modelfun(mean_time_wt,theta,[],[],[]);

        model_wt=N(:,1:2);


        [~,~,~,N] = modelfun(mean_time_ko,theta,[],[],[]);

        model_ko=N(:,3:4);



        %%

        subplot(1,4,1)
        hold on
        errorbar(mean_time_wt,measured(:,2),error(:,2),'ob','markerfacecolor','b')
        plot(mean_time_wt,model_wt(:,2),'r','linewidth',2)

        if model == best_model
            ylabel('best')

        end

        if model == 1
            title('spliced wt')
        elseif model == 8
            xlabel('pseudotime')
        end

        %
        subplot(1,4,2)
        hold on
        errorbar(mean_time_wt,measured(:,1),error(:,1),'ob','markerfacecolor','b')
        plot(mean_time_wt,model_wt(:,1),'r','linewidth',2)

        if model == 1
            title('unspliced wt')
        elseif model == 8
            xlabel('pseudotime')
        end



        %

        subplot(1,4,3)
        hold on
        errorbar(mean_time_ko,measured(:,4),error(:,4),'ob','markerfacecolor','b')
        plot(mean_time_ko,model_ko(:,2),'r','linewidth',2)



        if model == 1
            title('spliced ko')
        elseif model == 8
            xlabel('pseudotime')
        end

        %
        subplot(1,4,4)
        hold on
        errorbar(mean_time_ko,measured(:,3),error(:,3),'ob','markerfacecolor','b')
        plot(mean_time_ko,model_ko(:,1),'r','linewidth',2)


        if model == 1
            title('unspliced ko')
        elseif model == 8
            xlabel('pseudotime')
        end


        %
    end


    

    set(gcf, 'PaperUnits', 'centimeters');
    exportfig(gcf,['figure_model/' gene_name '_best.eps'],'FontMode', 'fixed','Fontsize',20,'color', 'cmyk','width',35,'height',8,'Renderer','painters','Lockaxes',0);%
end
end