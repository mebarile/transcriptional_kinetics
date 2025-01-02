function [negLogL,grad] = llKin_model6_all_samples(varargin)


if nargin >= 3
    theta = varargin{1};
    modelFun = varargin{2};
    D = varargin{3};

else
    error('Not enough inputs!');
end

data_input = D.pop.mean;

% default options
options.n_exp = 1;
if nargin == 4
    options = setdefault(varargin{4},options);
end

% options.alpha

options_sim.sensi = 1;
options_sim.maxsteps = 1e6;

% observation times
t = D.pop.t;

[status,~,~,u,~,us] = modelFun(t,theta,[],[],options_sim);

if status ~= 0
    error([],'AMICI simulation failed')
end
% check if solution is negative
if max(max(u<0))
    error('negative solution')
end

ind_w1(1,1) = 1;
ind_k1(1,1) = 1;

ind_w2(1,1) = 1;
ind_k2(1,1) = 1;

ind_w3(1,1) = 1;
ind_k3(1,1) = 1;


time_all = [D.pop.tw,D.pop.tm];

for iHelp = 2:20

    ind_w1(iHelp,1) = find(D.pop.t == time_all(iHelp,1));
    ind_k1(iHelp,1) = find(D.pop.t == time_all(iHelp,4));

    ind_w2(iHelp,1) = find(D.pop.t == time_all(iHelp,2));
    ind_k2(iHelp,1) = find(D.pop.t == time_all(iHelp,5));

    ind_w3(iHelp,1) = find(D.pop.t == time_all(iHelp,3));
    ind_k3(iHelp,1) = find(D.pop.t == time_all(iHelp,6));

end


N = [u(ind_w1,1:2),u(ind_w2,3:4),u(ind_w3,5:6),u(ind_k1,7:8),u(ind_k2,9:10),u(ind_k3,11:12)];

dNdtheta = [us(ind_w1,1:2,:),us(ind_w2,3:4,:),us(ind_w3,5:6,:),us(ind_k1,7:8,:),us(ind_k2,9:10,:),us(ind_k3,11:12,:)];

negLogL = 0;
grad = zeros(size(theta));


% negative log-likelihood for log pop size
for it = 1:20
    for col = 1:12

        err_idx1 = floor((col-1)/6);
        err_idx2 = 2-mod(col,2);

        err_idx = 2*err_idx1+err_idx2;


        erri = exp(theta(24+8+8+4*3+err_idx));

        negLogL = negLogL + 0.5*log(2*pi*erri^2/options.n_exp) + ...
            0.5*(data_input(it,col)-N(it,col)).^2/...
            (erri^2/options.n_exp);

        vec_add = zeros(40+12+4,1);

        vec_add(24+8+8+4*3+err_idx) = 1 - (data_input(it,col)-N(it,col))^2/erri^2;


        grad = grad  - ((data_input(it,col)-N(it,col))/...
            (erri^2/options.n_exp)*...
            squeeze(dNdtheta(it,col,:))')'+vec_add;


    end
end

%%regularization


negLogL = negLogL + options.alpha*(sum((theta(2:8)-theta(1:7)).^2)+...
    sum((theta(10:16)-theta(9:15)).^2)+ sum((theta(18:24)-theta(17:23)).^2)+...
    (sum((theta(26:32)-theta(25:31)).^2)+...
    sum((theta(34:40)-theta(33:39)).^2)));


grad = grad + options.alpha*2*([-(theta(2:8)-theta(1:7));0;-(theta(10:16)-theta(9:15));0;-(theta(18:24)-theta(17:23));0;
    -(theta(26:32)-theta(25:31));0;-(theta(34:40)-theta(33:39));0;zeros(16,1)]...
    +[0;theta(2:8)-theta(1:7);0;theta(10:16)-theta(9:15);0;theta(18:24)-theta(17:23);...
    0;theta(26:32)-theta(25:31);0;theta(34:40)-theta(33:39);zeros(16,1)]);

end
