function model = kinetics_model6_syms( )
% compartment model for the rna velocity equations

model.param = 'log';
model.forward = true;
model.adjoint = false;

%% STATES
% create state syms: two compartments
x = sym('x',[4,1]);

%% PARAMETERS ( for these sensitivities will be computed )

% create parameter syms
syms aw1 aw2 aw3 aw4 aw5 aw6 aw7 aw8 bw1 bw2 bw3 bw4 bw5 bw6 bw7 bw8 gw1 gw2 gw3 gw4 gw5 gw6 gw7 gw8 ak1 ak2 ak3 ak4 ak5 ak6 ak7 ak8 gk1 gk2 gk3 gk4 gk5 gk6 gk7 gk8 Iuw Isw Iuk Isk e1 e2 e3 e4

p = [aw1,aw2,aw3,aw4,aw5,aw6,aw7,aw8,bw1,bw2,bw3,bw4,bw5,bw6,bw7,bw8,gw1,gw2,gw3,gw4,gw5,gw6,gw7,gw8,ak1,ak2,ak3,ak4,ak5,ak6,ak7,ak8,gk1,gk2,gk3,gk4,gk5,gk6,gk7,gk8,Iuw,Isw,Iuk,Isk,e1,e2,e3,e4];

%% CONSTANTS ( for these no sensitivities will be computed )
% this part is optional and can be ommited

% create parameter syms
% k = sym('k',[2,1]);

%% SYSTEM EQUATIONS

% create symbolic variable for time
syms t
%
t_all =  [ 0    0.0682    0.1363    0.2045    0.2727    0.3409    0.4090    0.4772];

lt = length(t_all);



alpha_w = am_spline_try(t,lt,t_all(1),aw1,t_all(2),aw2,t_all(3),aw3,t_all(4),aw4,t_all(5),aw5,...
    t_all(6),aw6,t_all(7),aw7,t_all(8),aw8,0,0.0);

beta_w =  am_spline_try(t,lt,t_all(1),bw1,t_all(2),bw2,t_all(3),bw3,t_all(4),bw4,t_all(5),bw5,...
    t_all(6),bw6,t_all(7),bw7,t_all(8),bw8,0,0.0);

gamma_w = am_spline_try(t,lt,t_all(1),gw1,t_all(2),gw2,t_all(3),gw3,t_all(4),gw4,t_all(5),gw5,...
    t_all(6),gw6,t_all(7),gw7,t_all(8),gw8,0,0.0);

%

alpha_k = am_spline_try(t,lt,t_all(1),ak1,t_all(2),ak2,t_all(3),ak3,t_all(4),ak4,t_all(5),ak5,...
    t_all(6),ak6,t_all(7),ak7,t_all(8),ak8,0,0.0);
   


gamma_k = am_spline_try(t,lt,t_all(1),gk1,t_all(2),gk2,t_all(3),gk3,t_all(4),gk4,t_all(5),gk5,...
    t_all(6),gk6,t_all(7),gk7,t_all(8),gk8,0,0.0);


%
xdot = sym(zeros(size(x)));

xdot(1) = alpha_w - beta_w * x(1);
xdot(2) = beta_w * x(1) - gamma_w * x(2);
xdot(3) = alpha_k - beta_w * x(3);
xdot(4) = beta_w * x(3) - gamma_k * x(4);

%% INITIAL CONDITIONS



x0 = [Iuw,Isw,Iuk,Isk];

%% OBSERVALES
% N = sum(x);
% y = [x;N];

y = x;
%% SYSTEM STRUCT

model.sym.x = x;

model.sym.xdot = xdot;
model.sym.p = p;
model.sym.x0 = x0;
model.sym.y = y;
end
