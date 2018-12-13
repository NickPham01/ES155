clear all; clc; close all;

% Define parameters
sigma = 10;
beta = 8/3;
rho = 28;

A = [-10 10 0;
    8/3 -1 -1;
    0 0 28];

evalues = eig(A);

% Generate true state
N = 10000;
dt = 0.01;
x0 = [mvnrnd(0,1); mvnrnd(0,1); mvnrnd(0,1)];
xtrue = zeros(3,N);
xtrue(:,1) = x0;

for i = 2:N
    xtrue(1,i) = xtrue(1,i-1) + (sigma * (xtrue(2,i-1) - xtrue(1,i-1))) * dt;
    xtrue(2,i) = xtrue(2,i-1) + (rho * xtrue(1,i-1) - xtrue(2,i-1) - (xtrue(1,i-1) * xtrue(3,i-1))) * dt;
    xtrue(3,i) = xtrue(3,i-1) + (-beta * xtrue(3,i-1) + (xtrue(1,i-1) * xtrue(2,i-1))) * dt;
end

% Plot function
figure();
hold on;
plot3(xtrue(1,:), xtrue(2,:),xtrue(3,:),'LineWidth',1.2);
xlabel('x');
ylabel('y');
zlabel('z');
title('Lorenz System');

