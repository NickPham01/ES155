clear; clc;

d   = 1;
c   = 0.01;
w0  = 0.05;
l   = 0.5;

P = tf(d, [1, d]) * tf(c, [1, 2*l*w0, w0^2])
C = 20 * tf([11, 12, 1], [0.54, 1])
L = P*C
sys = feedback(L,1)

SSerr = 1 - dcgain(sys);
dbdrop = mag2db(0.9);
BW = bandwidth(sys, dbdrop);
[Gm,Pm,Wcg,Wcp] = margin(sys);

SSerr
BW
Pm
