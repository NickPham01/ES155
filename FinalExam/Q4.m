%% Problem 4

% part a)

% write the atomic transfer functions as
% u -> + -> M - - -> N -> y
%      |        |
%      --- -1 ---

a = 1

M = tf([1, a], [1, 19*a])
N = tf(1, [1, 0.02*a, 0.0001*a^2])

P = N * feedback(M, 1)

% check with
P_check = tf([1,a], [1, 10*a]) * tf(1, [1, 0.02*a, 0.0001]) * 1/2
isequal(P, P_check)

% part b)
figure(1); clf;

L = P*1
bode(L)
margin(L)