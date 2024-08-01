function write_all_parameters
%% This code writes the parameters as computed along a more dense pseudotime

vec_time  = 0:0.4772/7:0.4772;

%% read data

stat = dlmread(['./results/statistics_aic.txt']);

temp = [dlmread('./data_output/dpt_w.txt');dlmread('./data_output/dpt_m.txt')];
time_dpt = [temp;1:length(temp)];
time_dpt = sortrows(time_dpt',1)';
plot_time = time_dpt(1,:);
dlmwrite('./results/time_dpt_w_m.txt',time_dpt)


plot_time = plot_time*1.42;

%% write parameters

for n = 1:1228


    model = stat(n,2);

    all_parameters  = dlmread(['./results/fit_parameters.txt']);

    theta = all_parameters((n-1)*8+model,3:end);

    %% 

    alpha_w = spline(vec_time,theta(1:8),plot_time);
    beta_w = spline(vec_time,theta(9:16),plot_time);
    gamma_w = spline(vec_time,theta(17:24),plot_time);


    alpha_m = alpha_w;
    beta_m = beta_w;
    gamma_m = gamma_w;



    if model == 2

        alpha_m = spline(vec_time,theta(25:32),plot_time);

    elseif model == 3

        beta_m = spline(vec_time,theta(25:32),plot_time);

    elseif model == 4

        gamma_m = spline(vec_time,theta(25:32),plot_time);

    elseif model == 5

        alpha_m = spline(vec_time,theta(25:32),plot_time);
        beta_m = spline(vec_time,theta(33:40),plot_time);

    elseif model == 6

        alpha_m = spline(vec_time,theta(25:32),plot_time);
        gamma_m = spline(vec_time,theta(33:40),plot_time);


    elseif model == 7

        beta_m = spline(vec_time,theta(25:32),plot_time);
        gamma_m = spline(vec_time,theta(33:40),plot_time);

    elseif model == 8

        alpha_m = spline(vec_time,theta(25:32),plot_time);
        beta_m = spline(vec_time,theta(33:40),plot_time);
        gamma_m = spline(vec_time,theta(41:48),plot_time);

    end


    dlmwrite(['./results/all_parameters.txt'],[alpha_w-alpha_m,beta_w-beta_m,gamma_w-gamma_m],'-append')


end
end