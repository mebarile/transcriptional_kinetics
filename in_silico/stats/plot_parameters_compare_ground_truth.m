function plot_parameters_compare_ground_truth
%% this code compares the parameters estimated by model "nodel_fit" compared to the ground thruth


number_gene = 1;
model_ground_truth = 2;

model_fit = 1;

data = load(['./data_input2/data_gene_' num2str(number_gene) '_model' num2str(model_ground_truth) '.mat']);
%
para = dlmread('./reg10/parameters.txt');

line = find((para(:,1) == number_gene) & (para(:,2) == model_ground_truth) & (para(:,3) == model_fit));

theta_fit = para(line,4:end);



%%

vec_time  = 0:0.4772/7:0.4772;




dt = 0.0015;
plot_time = vec_time(1):dt:vec_time(end);

theta_real = data.D.pop.para_sim;
theta_real =theta_real';




%%


alfa_w_real = spline(vec_time,theta_real(1:8),plot_time);

beta_w_real = spline(vec_time,theta_real(9:16),plot_time);

gamma_w_real = spline(vec_time,theta_real(17:24),plot_time);



alfa_w_fit = spline(vec_time,theta_fit(1:8),plot_time);

beta_w_fit = spline(vec_time,theta_fit(9:16),plot_time);

gamma_w_fit = spline(vec_time,theta_fit(17:24),plot_time);
%%

if model_fit == 1


    alfa_m_fit = alfa_w_fit;

    beta_m_fit = beta_w_fit;

    gamma_m_fit = gamma_w_fit;


elseif model_fit == 2

    alfa_m_fit  = spline(vec_time,theta_fit(25:32),plot_time);

    beta_m_fit = beta_w_fit;

    gamma_m_fit = gamma_w_fit;


elseif model_fit == 3

    alfa_m_fit = alfa_w_fit;

    beta_m_fit = spline(vec_time,theta_fit(25:32),plot_time);

    gamma_m_fit = gamma_w_fit;


elseif model_fit == 4

    alfa_m_fit = alfa_w_fit;

    beta_m_fit = beta_w_fit;

    gamma_m_fit = spline(vec_time,theta_fit(25:32),plot_time);


elseif model_fit == 5

    alfa_m_fit  = spline(vec_time,theta_fit(25:32),plot_time);

    beta_m_fit = spline(vec_time,theta_fit(33:40),plot_time);

    gamma_m_fit = gamma_w_fit;

elseif model_fit == 6


    alfa_m_fit = spline(vec_time,theta_fit(25:32),plot_time);

    beta_m_fit = beta_w_fit;

    gamma_m_fit = spline(vec_time,theta_fit(33:40),plot_time);

elseif model_fit == 7


    alfa_m_fit = alfa_w_fit;

    beta_m_fit = spline(vec_time,theta_fit(25:32),plot_time);

    gamma_m_fit = spline(vec_time,theta_fit(33:40),plot_time);

elseif model_fit == 8

    alfa_m_fit  = spline(vec_time,theta_fit(25:32),plot_time);

    beta_m_fit = spline(vec_time,theta_fit(33:40),plot_time);

    gamma_m_fit = spline(vec_time,theta_fit(41:48),plot_time);


end


%%

if model_ground_truth == 1


    alfa_m_real = alfa_w_real;

    beta_m_real = beta_w_real;

    gamma_m_real = gamma_w_real;



elseif model_ground_truth == 2

    alfa_m_real  = spline(vec_time,theta_real(25:32),plot_time);

    beta_m_real = beta_w_real;

    gamma_m_real = gamma_w_real;


elseif model_ground_truth == 3

    alfa_m_real = alfa_w_real;

    beta_m_real = spline(vec_time,theta_real(25:32),plot_time);

    gamma_m_real = gamma_w_real;

elseif model_ground_truth == 4

    alfa_m_real = alfa_w_real;

    beta_m_real = beta_w_real;

    gamma_m_real = spline(vec_time,theta_real(25:32),plot_time);


elseif model_ground_truth == 5

    alfa_m_real  = spline(vec_time,theta_real(25:32),plot_time);

    beta_m_real = spline(vec_time,theta_real(33:40),plot_time);

    gamma_m_real = gamma_w_real;

elseif model_ground_truth == 6


    alfa_m_real = spline(vec_time,theta_real(25:32),plot_time);

    beta_m_real = beta_w_real;

    gamma_m_real = spline(vec_time,theta_real(33:40),plot_time);

elseif model_ground_truth == 7


    alfa_m_real = alfa_w_real;

    beta_m_real = spline(vec_time,theta_real(25:32),plot_time);

    gamma_m_real = spline(vec_time,theta_real(33:40),plot_time);

elseif model_ground_truth == 8

    alfa_m_real  = spline(vec_time,theta_real(25:32),plot_time);

    beta_m_real = spline(vec_time,theta_real(33:40),plot_time);

    gamma_m_real = spline(vec_time,theta_real(41:48),plot_time);
end







%%
figure(1)

clf

subplot(2,3,1)

hold on

plot(plot_time/max(plot_time),exp(alfa_w_fit),'linewidth',2,'color','r')
plot(plot_time/max(plot_time),exp(alfa_w_real),'linewidth',2,'color','b')
title('transcription')

ylabel('condition 1')
legend('fit','simulation')



subplot(2,3,2)
hold on
plot(plot_time/max(plot_time),exp(beta_w_fit),'linewidth',2,'color','r')
plot(plot_time/max(plot_time),exp(beta_w_real),'linewidth',2,'color','b')
title('splicing')


subplot(2,3,3)
hold on
plot(plot_time/max(plot_time),exp(gamma_w_fit),'linewidth',2,'color','r')
plot(plot_time/max(plot_time),exp(gamma_w_real),'linewidth',2,'color','b')
title('degradation')




subplot(2,3,4)

hold on

plot(plot_time/max(plot_time),exp(alfa_m_fit),'linewidth',2,'color','r')
plot(plot_time/max(plot_time),exp(alfa_m_real),'linewidth',2,'color','b')
title('transcription')


ylabel('condition 2')


subplot(2,3,5)
hold on
plot(plot_time/max(plot_time),exp(beta_m_fit),'linewidth',2,'color','r')
plot(plot_time/max(plot_time),exp(beta_m_real),'linewidth',2,'color','b')
xlabel('pseudotime')
title('splicing')



subplot(2,3,6)
hold on
plot(plot_time/max(plot_time),exp(gamma_m_fit),'linewidth',2,'color','r')
plot(plot_time/max(plot_time),exp(gamma_m_real),'linewidth',2,'color','b')
title('degradation')



end
