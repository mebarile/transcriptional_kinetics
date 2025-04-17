function plot_all_parameters_best
%% read data


genes = readtable('./data_input/list_genes.csv');


pop1 = 1;
pop2 = 2;

stat = dlmread(['./results/statistics_aic_' num2str(pop1) '_' num2str(pop2) '.txt']);

vec_time  = 0:0.4772/7:0.4772;
dt = 0.001;
plot_time = vec_time(1):dt:vec_time(end);

for n = 441 % gene number

    temp = genes.Gene(n);

    gene_name = temp{1}

    best_model = stat(n,2);

    all_para  = dlmread(['./results/statistics_aic_' num2str(pop1) '_' num2str(pop2) '.txt']);

    %%
    figure(11)
    clf


    for model = best_model

        theta = all_para((n-1)*8+model,3:end);

        %% model

        alfa_w = spline(vec_time,theta(1:8),plot_time);
        beta_w = spline(vec_time,theta(9:16),plot_time);
        gamma_w = spline(vec_time,theta(17:24),plot_time);


        alfa_m = alfa_w;
        beta_m = beta_w;
        gamma_m = gamma_w;



        if model == 2
            alfa_m = spline(vec_time,theta(25:32),plot_time);

        elseif model == 3


            beta_m = spline(vec_time,theta(25:32),plot_time);
        elseif model == 4

            gamma_m = spline(vec_time,theta(25:32),plot_time);
        elseif model == 5
            alfa_m = spline(vec_time,theta(25:32),plot_time);

            beta_m = spline(vec_time,theta(33:40),plot_time);


        elseif model == 6
            alfa_m = spline(vec_time,theta(25:32),plot_time);
            gamma_m = spline(vec_time,theta(33:40),plot_time);


        elseif model == 7

            beta_m = spline(vec_time,theta(25:32),plot_time);
            gamma_m = spline(vec_time,theta(33:40),plot_time);

        elseif model == 8

            alfa_m = spline(vec_time,theta(25:32),plot_time);
            beta_m = spline(vec_time,theta(33:40),plot_time);
            gamma_m = spline(vec_time,theta(41:48),plot_time);

        end


        %%


        subplot(1,3,1)
        hold on
        % plot(plot_time,alfa_w,'linewidth',2,'color','k')
        plot(plot_time,alfa_m - alfa_w,'linewidth',4,'color','b')
        line([0 plot_time(end)],[0 0],'color','k','linewidth',4)


        title('transcription')



        if model == best_model
            ylabel(['best model: ' num2str(best_model)])

        end


        %

        subplot(1,3,2)
        hold on
        plot(plot_time,beta_m - beta_w,'linewidth',4,'color','b')
        line([0 plot_time(end)],[0 0],'color','k','linewidth',4)


        title('splicing')
        xlabel('pseudotime')


        subplot(1,3,3)
        hold on
        plot(plot_time,gamma_m - gamma_w,'linewidth',4,'color','b')
        line([0 plot_time(end)],[0 0],'color','k','linewidth',4)


        title('degradation')

    end



    set(gcf, 'PaperUnits', 'centimeters');
    exportfig(gcf,['figures_parameters/' gene_name '_dif.eps'],'FontMode', 'fixed','Fontsize',20,'color', 'cmyk','width',35,'height',8,'Renderer','painters','Lockaxes',0);%
end
end