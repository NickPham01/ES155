%% Problem 1

% Part c

% parameters (arbitrary satisfy given assumptions ( > 0)
%syms omega_0 zeta_0 omega_e zeta_e


omega_0 = 0.1
zeta_0 = 0.2
omega_e = 0.3
zeta_e = 0.4

% constants
k1 = 4*zeta_0*omega_0 - 8
k2 = 2*omega_0^2 - 4*zeta_0*omega_0 + 6
kr = 2*omega_0^2

l1 = omega_e^2 - 6*zeta_e*omega_e + 23
l2 = 2*zeta_e*omega_e - 4

% Matrices
A = [-3 2; 1 -1]
B = [0.5; 0]
C = [1 0]

K = [k1 k2]
L = [l1; l2]

% Compute the eigenvalues to 
%A_tilde = [(A - B*K), zeros(2); zeros(2), (A - L*C)]
A_tilde = [(A - B*K), B*K; zeros(2), (A - L*C)]
eig(A_tilde)
%latex(eig(A_tilde))
