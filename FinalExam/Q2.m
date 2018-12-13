%% Problem 2
clear; clc;

%% 2.a
% Give constants arbitrary values
syms C0 C1 a b

% Write the state space matrices
A = [-(a+b)/C0, a/C0; a/C1, -a/C1]
B = [1/C0; 0]

% Check if the system is stable
eig(A)

C = [1, 0; 0, 1]
D = [0; 0]


%% 2.b

Wr = [B A*B]
rank(Wr)

syms k1 k2
K = [k1 k2]
eig(A - B*K)

syms C

place(A, B, [-b/C, -b/C])