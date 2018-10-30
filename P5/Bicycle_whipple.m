% bicycle_whipple.m
% KJA, 20 Aug 07
%
% This file contains the parameters that are used for the Whipple
% bicycle model, introduced in Section 3.2 of AM08.  The model is
% based on the linearized 4th order model and analysis of eigenvalues
% from IEEE CSM (25:4) August 2005 pp 26-47

clear;

%% Given:
%   Basic data is given by 26 parameters
g = 9.81;			% Acceleration of gravity [m/s^2]
b = 1.00;			% Wheel base [m]
c = 0.08;			% Trail [m]
Rrw = 0.35; Rfw = 0.35;		% Wheel radii
lambda = pi*70/180;		% Head angle [radians]

% Rear frame mass [kg], center of mass [m], and inertia tensor [kgm^2]
mrf=12;xrf=0.439;zrf=0.579;
Jxxrf=0.475656;Jxzrf=0.273996;Jyyrf=1.033092;Jzzrf=0.527436;
mrf=87;xrf=0.491586;zrf=1.028138;
Jxxrf=3.283666;Jxzrf=0.602765;Jyyrf=3.8795952;Jzzrf=0.565929;

% Front frame mass [kg], center of mass [m], and inertia tensor [kgm^2]
mff=2;xff=0.866;zff=0.676;
Jxxff=0.08;Jxzff=-0.02;Jyyff=0.07;Jzzff=0.02;

% Rear wheel mass [kg], center of mass [m], and inertia tensor [kgm^2]
mrw=1.5;Jxxrw=0.07;Jyyrw=0.14;

% Front wheel mass [kg], center of mass [m], and inertia tensor [kgm^2]
mfw=1.5;Jxxfw=0.07;Jyyfw=0.14;

% Auxiliary variables
xrw=0;zrw=Rrw;xfw=b;zfw=Rfw;
Jzzrw=Jxxrw;Jzzfw=Jxxfw;
mt=mrf+mrw+mff+mfw;
xt=(mrf*xrf+mrw*xrw+mff*xff+mfw*xfw)/mt;
zt=(mrf*zrf+mrw*zrw+mff*zff+mfw*zfw)/mt;
Jxxt=Jxxrf+mrf*zrf^2+Jxxrw+mrw*zrw^2+Jxxff+mff*zff^2+Jxxfw+mfw*zfw^2;
Jxzt=Jxzrf+mrf*xrf*zrf+mrw*xrw*zrw+Jxzff+mff*xff*zff+mfw*xfw*zfw;
Jzzt=Jzzrf+mrf*xrf^2+Jzzrw+mrw*xrw^2+Jzzff+mff*xff^2+Jzzfw+mfw*xfw^2;
mf=mff+mfw;
xf=(mff*xff+mfw*xfw)/mf;zf=(mff*zff+mfw*zfw)/mf;
Jxxf=Jxxff+mff*(zff-zf)^2+Jxxfw+mfw*(zfw-zf)^2;
Jxzf=Jxzff+mff*(xff-xf)*(zff-zf)+mfw*(xfw-xf)*(zfw-zf);
Jzzf=Jzzff+mff*(xff-xf)^2+Jzzfw+mfw*(xfw-xf)^2;
d=(xf-b-c)*sin(lambda)+zf*cos(lambda);
Fll=mf*d^2+Jxxf*cos(lambda)^2+2*Jxzf*sin(lambda)*cos(lambda)+Jzzf*sin(lambda)^2;
Flx=mf*d*zf+Jxxf*cos(lambda)+Jxzf*sin(lambda);
Flz=mf*d*xf+Jxzf*cos(lambda)+Jzzf*sin(lambda);
gamma=c*sin(lambda)/b;
Sr=Jyyrw/Rrw;Sf=Jyyfw/Rfw;St=Sr+Sf;Su=mf*d+gamma*mt*xt;

% Matrices for linearized fourth order model
M=[Jxxt -Flx-gamma*Jxzt;-Flx-gamma*Jxzt Fll+2*gamma*Flz+gamma^2*Jzzt];
K0=[-mt*g*zt g*Su;g*Su  -g*Su*cos(lambda)];
K2=[0 -(St+mt*zt)*sin(lambda)/b;0 (Su+Sf*cos(lambda))*sin(lambda)/b];
c12=gamma*St+Sf*sin(lambda)+Jxzt*sin(lambda)/b+gamma*mt*zt;
c22=Flz*sin(lambda)/b+gamma*(Su+Jzzt*sin(lambda)/b);
C0=[0 -c12;(gamma*St+Sf*sin(lambda)) c22]; 
one=diag([1 1]);null=zeros(2,2);

% Nominal velocity 
v0=5;

% Matrices of state model
A=[null one;-M\(K0+K2*v0^2) -M\(C0*v0)];
bm=M\[0;1];
B=[0;0;bm];    
eig(A)';


%% Compute K for different given eigenvalues

%{
eigs = [-2, -10, -1+i, -1-i;
        -2, -10, -2+2i, -2-2i;
        -2, -10, -5+5i, -5-5i];
    
C = [0 1 0 0];
D = 0;

figure(2); clf; hold on;
for i = 1:size(eigs,1)
    fprintf("For the eigenvalues:\n")
    eigs(i,:)
    
   K = place(A, B, eigs(i,:))
   kr = inv(-(C - D*K)*inv(A - B*K)*B + D)
   
   sys = ss(A-B*K, kr*B, C, D);
   opt = stepDataOptions('StepAmplitude',0.002);
   [y, t, x] = step(sys, 6, opt);
   figure(2)
   subplot(2,1,1)
   hold on;
   plot(t, y)
   ylim([-1, 2.5] * 10^-3)
   
   subplot(2,1,2)
   hold on;
   T = -K*x' + kr*0.002;
   plot(t, T)
   ylim([-0.03, 0.005])
   
end
hold off

eigs_legend = ({"$\lambda = -2, -10, -1 \pm i$", "$\lambda = -2, -10, -2 \pm 2i$", "$\lambda = -2, -10, -5 \pm 5i$"});
titles = {"Output Steering Angle $\delta$, radians", "Input Torque T"};
for i = 1:2
    subplot(2,1,i);
    legend(eigs_legend, 'Interpreter', 'latex', 'Location', 'southeast')
    xlabel("Time $(s)$", 'Interpreter', 'latex')
    title(titles(i), 'Interpreter', 'latex')
end

saveas(gca, "ES155P4_2_bicycleStepResponse.jpg")

%}

%% Homework 5.2
%% 2.a
fprintf(['Part 2.a\n'])
A
B
C = [1 0 0 0]

w_o = [C; C*A; C*A^2; C*A^3]
rank(w_o)

%% 2.b
fprintf(['Part 2.b\n'])

L_eigs = [-4, -20, -2 + 2i, -2 - 2i];
LT = place(A', C', L_eigs);
L = LT'

%% 2.c
fprintf(['Part 2.c\n'])

D = 0
K_eigs = [-2, -10, -1+i, -1 - i]
K = place(A, B, K_eigs)
kr = inv(-(C - D*K)*inv(A - B*K)*B + D)

A_t = [A-B*K, B*K; zeros(size(A)), A - L*C]
B_t = [B*kr; zeros(size(B))]
C_t = [C zeros(size(C))]

% simulate
sys = ss(A_t, B_t, C_t, 0);
opt = stepDataOptions('StepAmplitude',0.002);
[y, t, x] = step(sys, 6, opt);

% plot the output tilt angle
figure(1); clf;
subplot(2,1,1);
plot(t, y);
title('Tilt Angle $y = \phi$', 'Interpreter', 'latex')
xlabel('$t$', 'Interpreter', 'latex')
ylabel('$y = \phi$', 'Interpreter', 'latex')

% compute x_hat by extracting the state x and error e
e = x(:, 5:8);
x = x(:, 1:4);
x_hat = x + e;
% 
% get the Torque input then plot
T = -K*x_hat' + kr*0.002;
subplot(2,1,2);
plot(t, x(:,2));
title('Torque $y = \phi$', 'Interpreter', 'latex')
xlabel('$t$', 'Interpreter', 'latex')
ylabel('$y = \phi$', 'Interpreter', 'latex')

saveas(gca, "ES155P5_2_bicycle.jpg")


