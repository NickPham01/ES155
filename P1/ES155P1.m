%% Problem 4: Cruise Control

%% Part A: Equilibrium for open and closed loop controllers

% beta = a/a_hat: use as domain from 0:2
beta = 0:0.05:2

% let Vss/Vref = V
V_ff = 1./beta
V_fb = 10./(beta + 10)

figure(1)
clf
hold on
plot(beta, V_ff)
plot(beta, V_fb)

legend(["V_ff", "V_fb"])
title("Steady State Velocity vs \beta")
xlabel("\beta")
ylabel(['$$v$$ steady state'], 'interpreter', 'latex')
hold off