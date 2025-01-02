function model = kinetics_model1_syms( )
% compartment model for the rna velocity equations

model.param = 'log';
model.forward = true;
model.adjoint = false;

%% STATES
x = sym('x',[12,1]);

%% PARAMETERS ( for these sensitivities will be computed )

% create parameter syms

syms aw1 aw2 aw3 aw4 aw5 aw6 aw7 aw8 bw1 bw2 bw3 bw4 bw5 bw6 bw7 bw8 gw1 gw2 gw3 gw4 gw5 gw6 gw7 gw8 ak1 ak2 ak3 ak4 ak5 ak6 ak7 ak8 bk1 bk2 bk3 bk4 bk5 bk6 bk7 bk8 Iuw1 Isw1 Iuw2 Isw2 Iuw3 Isw3 Iuk1 Isk1 Iuk2 Isk2 Iuk3 Isk3 e1 e2 e3 e4

p = [aw1,aw2,aw3,aw4,aw5,aw6,aw7,aw8,bw1,bw2,bw3,bw4,bw5,bw6,bw7,bw8,gw1,gw2,gw3,gw4,gw5,gw6,gw7,gw8,ak1,ak2,ak3,ak4,ak5,ak6,ak7,ak8,bk1,bk2,bk3,bk4,bk5,bk6,bk7,bk8,Iuw1,Isw1,Iuw2,Isw2,Iuw3,Isw3,Iuk1,Isk1,Iuk2,Isk2,Iuk3,Isk3,e1,e2,e3,e4];



%% SYSTEM EQUATIONS

% create symbolic variable for time
syms t
% 
t_all =  [ 0    0.0682    0.1363    0.2045    0.2727    0.3409    0.4090    0.4772]; % adjust this on your dpt length
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

beta_k =  am_spline_try(t,lt,t_all(1),bk1,t_all(2),bk2,t_all(3),bk3,t_all(4),bk4,t_all(5),bk5,...
    t_all(6),bk6,t_all(7),bk7,t_all(8),bk8,0,0.0);

%
xdot = sym(zeros(size(x)));

xdot(1) = alpha_w - beta_w * x(1);
xdot(2) = beta_w * x(1) - gamma_w * x(2);
xdot(3) = alpha_w - beta_w * x(3);
xdot(4) = beta_w * x(3) - gamma_w * x(4);
xdot(5) = alpha_w - beta_w * x(5);
xdot(6) = beta_w * x(5) - gamma_w * x(6);
xdot(7) = alpha_k - beta_k * x(7);
xdot(8) = beta_k * x(7) - gamma_w * x(8);
xdot(9) = alpha_k - beta_k * x(9);
xdot(10) = beta_k * x(9) - gamma_w * x(10);
xdot(11) = alpha_k - beta_k * x(11);
xdot(12) = beta_k * x(11) - gamma_w * x(12);

%% INITIAL CONDITIONS

x0 = [Iuw1,Isw1,Iuw2,Isw2,Iuw3,Isw3,Iuk1,Isk1,Iuk2,Isk2,Iuk3,Isk3];

%% OBSERVALES


y = x;
%% SYSTEM STRUCT

model.sym.x = x;

model.sym.xdot = xdot;
model.sym.p = p;
model.sym.x0 = x0;
model.sym.y = y;
end

