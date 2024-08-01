function [negLogL,grad] = llKin_model8(varargin)
% compute the log likelihood value for Qiu compartment data for constant
% transition and time-dependent growth rates




if nargin >= 3
    theta = varargin{1};
    modelFun = varargin{2};
    D = varargin{3};
    pop1 = varargin{4};
    pop2 = varargin{5};
    
    
else
    error('Not enough inputs!');
end

data_input = D.pop.mean;
data_input = data_input(:,[(pop1-1)*2+1,pop1*2,(pop2-1)*2+1,pop2*2]);


% default options
options.n_exp = 1;
if nargin == 6
    options = setdefault(varargin{6},options);
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


ind_w(1,1) = 1;
ind_k(1,1) = 1;

time_all = [D.pop.tw,D.pop.tm];

for iHelp = 2:20
    
    ind_w(iHelp,1) = find(D.pop.t == time_all(iHelp,pop1));
    ind_k(iHelp,1) = find(D.pop.t == time_all(iHelp,pop2));
    
end


N = [u(ind_w,1:2),u(ind_k,3:4)];

% size(N)
dNdtheta = [us(ind_w,1:2,:),us(ind_k,3:4,:)];
% size(dNdtheta)


negLogL = 0;
grad = zeros(size(theta));


% negative log-likelihood for log pop size
for it = 1:20
    for col = 1:4
        
        erri = exp(theta(52+col));
        
        negLogL = negLogL + 0.5*log(2*pi*erri^2/options.n_exp) + ...
            0.5*(data_input(it,col)-N(it,col)).^2/...
            (erri^2/options.n_exp);
        
        vec_add = zeros(56,1);
        
        vec_add(52+col) = 1 - (data_input(it,col)-N(it,col))^2/erri^2;
        
        grad = grad  - ((data_input(it,col)-N(it,col))/...
            (erri^2/options.n_exp)*...
            squeeze(dNdtheta(it,col,:))')'+vec_add;
        
    end
end

%%regularization


negLogL = negLogL + options.alpha*(sum((theta(2:8)-theta(1:7)).^2)+...
    sum((theta(10:16)-theta(9:15)).^2)+ sum((theta(18:24)-theta(17:23)).^2)+...
    (sum((theta(26:32)-theta(25:31)).^2)+...
    sum((theta(34:40)-theta(33:39)).^2)+ sum((theta(42:48)-theta(41:47)).^2)));



grad = grad + options.alpha*2*([-(theta(2:8)-theta(1:7));0;-(theta(10:16)-theta(9:15));0;-(theta(18:24)-theta(17:23));0;
                                -(theta(26:32)-theta(25:31));0;-(theta(34:40)-theta(33:39));0;-(theta(42:48)-theta(41:47));0;0;0;0;0;0;0;0;0]...
    +[0;theta(2:8)-theta(1:7);0;theta(10:16)-theta(9:15);0;theta(18:24)-theta(17:23);...
    0;theta(26:32)-theta(25:31);0;theta(34:40)-theta(33:39);0;theta(42:48)-theta(41:47);0;0;0;0;0;0;0;0]);
%

end
