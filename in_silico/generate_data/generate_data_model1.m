function generate_data_model1
%% This code generates the genes whose kinetic rates follow Model 1 rules.
% For each of the 100 genes, the code creates a data structure D, which contains the following parts:
% D.pop.tw: the mean pseudotime of the cells in condition w (w is usually the wild type)
% D.pop.tm: the mean pseudotime of the cells in condition m (m is usually the mutant, KO, etc.)

% D.pop.t: the overall mean pseudotime, sorted
% D.pop.mean: mean spliced and unspliced counts in both conditions
% D.pop.var: variance of the spliced and unspliced counts in both conditions
% D.pop.para_sim: the kinetic rates used to produced the counts, and the
% error on the rates, also generated


%% General definitions
scale_time = 1;
interval = [1.1;0.9];

width = 0.2;
n_iterations = 20000;

nbins = 20;
vec_time  = 0:0.4772/7:0.4772;

mean_time = 0:0.001:vec_time(end);

%% Generate values

for gene = 1:100

    theta = zeros(32,1);
    
    sign_temp = rand(3,1);
    sign_temp(sign_temp >= 0.5) = 2;
    sign_temp(sign_temp < 0.5) = 1;

    %% generate the kinetic rates

    theta(1:8:24) = rand(3,1);
    theta(2:8:24) = randn(3,1) .* abs(theta(1:8:24)) * width + theta(1:8:24) .* interval(sign_temp);
    theta(3:8:24) = randn(3,1) .* abs(theta(2:8:24)) * width + theta(2:8:24) .* interval(sign_temp);
    theta(4:8:24) = randn(3,1) .* abs(theta(3:8:24)) * width + theta(3:8:24) .* interval(sign_temp);
    theta(5:8:24) = randn(3,1) .* abs(theta(4:8:24)) * width + theta(4:8:24) .* interval(sign_temp);
    theta(6:8:24) = randn(3,1) .* abs(theta(5:8:24)) * width + theta(5:8:24) .* interval(sign_temp);
    theta(7:8:24) = randn(3,1) .* abs(theta(6:8:24)) * width + theta(6:8:24) .* interval(sign_temp);
    theta(8:8:24) = randn(3,1) .* abs(theta(7:8:24)) * width + theta(7:8:24) .* interval(sign_temp);

    %% generate the error on the counts

    theta(25) = log(.1+3*rand);
    theta(26) = log(2+2*rand);
    theta(27) = log(.1+3*rand);
    theta(28) = log(2+2*rand);


    %% generate the spliced and unspliced counts

    [~,~,~,N,~,~] = simulate_kinetics_model1(mean_time,theta,[],[],[]);


    model_wt = N(:,1:2);
    model_ko = N(:,3:4);


    for i=1:n_iterations

        time_ind = ceil(rand * length(mean_time));

        model_ko_u_sim(i,:) = [time_ind,randn * model_ko(time_ind,1) * width+ model_ko(time_ind,1)];
        model_ko_s_sim(i,:) = [time_ind,randn * model_ko(time_ind,2) * width+ model_ko(time_ind,2)];




        model_wt_u_sim(i,:) = [time_ind,randn * model_wt(time_ind,1) * width+ model_wt(time_ind,1)];
        model_wt_s_sim(i,:) = [time_ind,randn * model_wt(time_ind,2) * width+ model_wt(time_ind,2)];

    end



    model_ko_u_sim = sortrows(model_ko_u_sim,1);
    model_ko_s_sim = sortrows(model_ko_s_sim,1);
    model_wt_u_sim = sortrows(model_wt_u_sim,1);
    model_wt_s_sim = sortrows(model_wt_s_sim,1);


    time_sim_w = mean_time(model_ko_u_sim(:,1));
    time_sim_m = mean_time(model_wt_u_sim(:,1));


    out_name =  ['./data_input/data_gene_' num2str(gene) '_model1.mat'];


    measured_w = [model_wt_u_sim(:,2),model_wt_s_sim(:,2),time_sim_w'];
    measured_m = [model_ko_u_sim(:,2),model_ko_s_sim(:,2),time_sim_m'];


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
    D.pop.para_sim = theta;

    %% saving
    save(out_name,'D');

end

end