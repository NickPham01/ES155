clear all; clc; close all;

% Define constants as given in the problem
b1 = .1;
J1 = 10;
r1 = 1;
b2 = .1;
J2 = 5;
r2 = 1;
kequ = 1;

% 2a - Define system matrices
A = [-b1/J1 0 r1/J1;
    0 -b2/J2 -r2/J2;
    -kequ*r1 -kequ*r2 0];

B = [1/J1 0;
    0 1/J2;
    0 0];

C = [0 1 0;
     1 0 0];

D = 0;

% 2b - Calculate open-loop eigenvalues
eig_ol = eig(A);

% 2c - Calculate controllability matrix and its rank
cont = ctrb(A, B);
cont_rank = rank(cont);

% 2d - Calculate observability matrix and its rank
obs = obsv(A, C);
obsv_rank = rank(obs);

% 2e - Design feedback controller
G = place(A,B,[-3,-4,-5]);
A_cl = A - B*G;
eig_cl = eig(A_cl);

% 2f - Design observer
C_new = [1 0 0];
L = place(A',C_new',[-10,-8,-6]);
A_err = A-L'*C_new;
% 2g - Simulate combined system
ti = 0;
tf = 10;
x0 = [0.5; 0.5; 1; mvnrnd(0,1,3)];

% Simulate
[tout, yout] = ode45(@(t,y)sim(t,y,A_cl,B,G,A_err),[ti tf], x0);

% Plot simulation results
figure();
subplot(3,1,1)
plot(tout,yout(:,1),'LineWidth',2);
ylabel('w_1');
subplot(3,1,2)
plot(tout,yout(:,2),'LineWidth',2);
ylabel('w_2');
subplot(3,1,3)
plot(tout,yout(:,3),'LineWidth',2);
ylabel('T');


function dydt = sim(~,y,A_cl,B,G,A_err)
% Use an augmented state and simulate the state+observer system
dydt = [A_cl -B*G;
        zeros(3,3) A_err]*y;
end