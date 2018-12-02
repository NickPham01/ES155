%% Given
s = tf('s')
m = 1, c=0.2, k=1
P = tf(1, [m,c,k])

figure(1);clf;
bode(P)
margin(P)

%% Design Control

% To acheive steady state error SSerr = 1/(1+L(0)) <= 0.01, want L(0) > 100
% 10% tracking error up to 1 rad/s so gain > 10 from 0-1 rad/s
% 45 degree phase margin

%controlSystemDesigner(P)

% From the controlSystemDesigner (pole at 1.65)
% C = tf(140) * tf([1, 0.6], 1)
%   The SSerr was just shy of reaching the goal (0.0118)
%   so try
%C = 200 * tf([1, 0.6], 1)
C = 150*tf([1, 1], 1)


sys = feedback(P*C,1)

figure(2); clf;
bode(sys)

[Gm, Pm, wgs, wps] = margin(sys);
Pm

SSerr = 1 - dcgain(sys)

