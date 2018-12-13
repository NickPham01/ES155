clear all; clc; close all;
set(0,'defaultLineLineWidth',2)

s = tf('s');

% System 1
G1 = (s+4)/(s*(s+10));
figure();
bode(G1)
figure();
step(G1/(1+G1)) %Closed loop step response
figure();
nyquist(G1)

% System 2
G2 = s/(s^2 + s + 4);
figure();
bode(G2);
figure();
step(G2/(1+G2))
figure();
nyquist(G2)

% System 3
G3 = (s+100)*s/((s+1000)*(s^2 - 4*s + 25));
figure();
bode(G3)
figure();
step(G3/(1+G3));
figure();
nyquist(G3)

set(findall(gcf,'type','line'),'linewidth',2)