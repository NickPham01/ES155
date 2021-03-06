%% Problem 1

%% First Examine System 1:
%% 1.1c

% let
g = 10

% Write the system
A = [0 1; -2*g 1]

eig(A)
% eigenvalues have positive reals so are not stable.  Check this with state
% space model:
B = [0;0]
C = [0 0]
D = 0

sys = ss(A,B,C,D)

isstable(sys)

%% Examine System 2:

A = [1 -1; 2 -1]
eig(A)

sys = ss(A,B,C,D)
isstable(sys)