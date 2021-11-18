clear all;
clc;




%%% Parameters Of The Quadcopter Vehicle:

m = 0.4080;                 % Kg
g = 9.8066;                 % m/s^2
d = 0.0301;                 % N.m/N
l = 0.1785;                 % m
Jx = 0.0022;                % Kg.m^2
Jy = 0.0029;                % Kg.m^2
Jz = 0.0048;                % Kg.m^2
fvmax = 1.93;               % N
Jrz = 2.029 * 10^(-5);      % Kg.m^2
Delta_q = [1 1 1 1;1 -1 -1 1;1 1 -1 -1;-1 1 -1 1];
Omega_M = 475;              % rad/s      %% The overall residual angular speed

%%% Parameters Of The PVTOL Vehicle:

Delta_p = [1 1;1 -1];

%%% Linear Dynamic model of the Quadcopter Around Hover Condition for ACAI:

Gc = [m*g 0 0 0]';
Jc = diag([1/m,l/Jx,l/Jy,d/Jz]);

Ac = [zeros(4,4) eye(4);zeros(4,4) zeros(4,4)];
Bc = [zeros(4,4);Jc];

%%% VTOL Controllability Analysis:

LoE = [0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9]';
ACAI_quad = [1.47 1.08 0.70 0.31 1.99 1.26 0.72 0.44 0.17]';
ACAI_PVTOL = [1.47 1.08 0.70 0.31 -0.06 -0.45 -0.84 -1.22 -1.61]';
LoE_ACAIq = [LoE,ACAI_quad];
LoE_ACAIp = [LoE,ACAI_PVTOL];
LoE_ACAIq_soft = LoE_ACAIq(1:4,:);
LoE_ACAIp_soft = LoE_ACAIp(1:4,:);
LoE_ACAIq_aggr = LoE_ACAIq(5:end,:);
LoE_ACAIp_aggr = LoE_ACAIp(5:end,:);

%%% FE using Nonlinear AO for PVTOL:

Ly1 = 0.009;
Ly2 = 5;
Lf1 = 50;
Lf2 = 50;

%%% FE using linear PIO for PVTOL:

Jp = diag([1/m,l/Jx]);
g_bold = [0,0,-g,0]';

A_breve = [zeros(2) eye(2);zeros(2) zeros(2)];
B_breve = [zeros(2);Jp];
C_breve = eye(4);
E_breve = [1 0;0 1;1 0;0 1];            % E has to be chosen
W_breve = [0 0 1 1]';                   % W has to be chosen

A_brevebar = [A_breve E_breve;zeros(2,4) zeros(2,2)];
B_brevebar = [B_breve;zeros(2,2)];
Gamma_brevebar = [W_breve;zeros(2,1)];
g_bar = [g_bold;zeros(2,1)];
C_brevebar = [C_breve;zeros(2,4)]';

% Finding the observer gain:

gamma = 0.064;

cvx_begin sdp


variable P(6,6) symmetric
variable M(6,4)
subject to

    [P*A_brevebar-M*C_brevebar+(P*A_brevebar-M*C_brevebar)'+eye(6) P*Gamma_brevebar;Gamma_brevebar'*P (-gamma^2)*eye(1)] <= 0; 

    P >= 0;

cvx_end

cvx_status

LPI_brevebar = (P^-1)*M;

%%% FE using qLPV PIO for Quadcopter:

Jr = diag([l/Jx,l/Jy,d/Jz]);
h1 = (Jy-Jz)/Jx;
h2 = (Jz-Jx)/Jy;
h3 = (Jx-Jy)/Jz;
b1 = -Jrz/Jx;
b2 = Jrz/Jy;
B = [zeros(3,3);Jr];
C = eye(6);
E = [1 0 0;0 1 0;0 0 1;1 0 0;0 1 0;0 0 1];        % E has to be chosen

kesi1 = -pi/2; kesi2 = -pi/2;
Ar = [0 0 h1*(kesi2);0 0 h2*(kesi1);h3*(-kesi2) 0 0];
A1 = [zeros(3) eye(3);zeros(3) Ar];
W1 = [zeros(1,3) b1*kesi2 b2*kesi1 0]';
A1_bar = [A1 E;zeros(3,6) zeros(3,3)];
Gamma1_bar = [W1;zeros(3,1)];

kesi1 = -pi/2; kesi2 = pi/2;
Ar = [0 0 h1*(kesi2);0 0 h2*(kesi1);h3*(-kesi2) 0 0];
A2 = [zeros(3) eye(3);zeros(3) Ar];
W2 = [zeros(1,3) b1*kesi2 b2*kesi1 0]';
A2_bar = [A2 E;zeros(3,6) zeros(3,3)];
Gamma2_bar = [W2;zeros(3,1)];

kesi1 = pi/2; kesi2 = -pi/2;
Ar = [0 0 h1*(kesi2);0 0 h2*(kesi1);h3*(-kesi2) 0 0];
A3 = [zeros(3) eye(3);zeros(3) Ar];
W3 = [zeros(1,3) b1*kesi2 b2*kesi1 0]';
A3_bar = [A3 E;zeros(3,6) zeros(3,3)];
Gamma3_bar = [W3;zeros(3,1)];

kesi1 = pi/2; kesi2 = pi/2;
Ar = [0 0 h1*(kesi2);0 0 h2*(kesi1);h3*(-kesi2) 0 0];
A4 = [zeros(3) eye(3);zeros(3) Ar];
W4 = [zeros(1,3) b1*kesi2 b2*kesi1 0]';
A4_bar = [A4 E;zeros(3,6) zeros(3,3)];
Gamma4_bar = [W4;zeros(3,1)];

B_bar = [B;zeros(3,3)];
C_bar = [C;zeros(3,6)]';

% Finding the observer gains:

gamma = 0.064;

% index 1

cvx_begin sdp

variable P(9,9) symmetric
variable M(9,6)

[P*A1_bar-M*C_bar+(P*A1_bar-M*C_bar)'+eye(9) P*Gamma1_bar;Gamma1_bar'*P (-gamma^2)*eye(1)] <= 0; 

P >= 0;

cvx_end

cvx_status

LPI1_bar = (P^-1)*M;

% index 2

cvx_begin sdp

variable P(9,9) symmetric
variable M(9,6)

[P*A2_bar-M*C_bar+(P*A2_bar-M*C_bar)'+eye(9) P*Gamma2_bar;Gamma2_bar'*P (-gamma^2)*eye(1)] <= 0; 

P >= 0;

cvx_end

cvx_status

LPI2_bar = (P^-1)*M;

% index 3

cvx_begin sdp

variable P(9,9) symmetric
variable M(9,6)

[P*A3_bar-M*C_bar+(P*A3_bar-M*C_bar)'+eye(9) P*Gamma3_bar;Gamma3_bar'*P (-gamma^2)*eye(1)] <= 0; 

P >= 0;

cvx_end

cvx_status

LPI3_bar = (P^-1)*M;

% index 4

cvx_begin sdp

variable P(9,9) symmetric
variable M(9,6)

[P*A4_bar-M*C_bar+(P*A4_bar-M*C_bar)'+eye(9) P*Gamma4_bar;Gamma4_bar'*P (-gamma^2)*eye(1)] <= 0; 

P >= 0;

cvx_end

cvx_status

LPI4_bar = (P^-1)*M;


%%% Nominal Controller:

K_P_kesi = diag([0.18,0.18,052]);
K_P_omega = diag([1.5,1.5,2]);
K_D_kesi = diag([0.08,0.08,0.22]);
K_D_omega = diag([0.2,0.2,0.5]);

Omega_d = [0 0 0]';
kesi_d = [0 0 0]';
O_d = [0 0 0]';
kesi_dot_d = [0 0 0]';

q_d = quaternion(Omega_d','rotvec');
q_d_conj = conj(q_d);


%%% Initial Conditions:

Phi_dot_init = 0;
Theta_dot_init = 0;
Psi_dot_init = 0;
x_dot_init = 0;
y_dot_init = 0;
z_dot_init = 0;

Phi_init = 0;
Theta_init = 0;
Psi_init = 0;
x_init = 0;
y_init = 0;
z_init = 0;


