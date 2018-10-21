%% ES155 Problem Set 4
%% 1.c

% define constants/parameters
w0 = 1;
z_0 = [0.1 0.4 0.7 0.9];
a0 = 1;
a1 = 2;
a2 = 1;
b0 = 0.5;

C = [0 1]
D = 0

eigs = zeros(2,length(z_0));

figure(1); clf;
hold on;

% try for each z0
for i = 1:length(z_0)
    z0 = z_0(i)
    
    % calculate control values
    k1 = 4*z0*w0 - 8;
    k2 = 2*w0^2 - 4*z0*w0 + 6;
    kr = 2*w0^2;

    % compute matrices
    A = [-a0 - a1, a1; a2, -a2];
    B = [b0; 0];
    K = [k1 k2];

    % compute the eigenvalues of (A - BK)
    eigs(:,i) = eig(A - B*K);
   
    % make state space model and plot step response
    sys = ss(A-B*K, kr*B, C, D);
    step(sys)
end
hold off;

z0_legend = strcat("${\zeta}_{0} = ", strtrim(cellstr(num2str(z_0'))'), "$")
legend(z0_legend, 'Interpreter', 'latex')

w = warning ('off','all');

fprintf(['The eigenvalues of the closed loop system response $(A - B*K)$ are '])
for i = 1:length(z_0)
    fprintf(['$%1.4f$ and $%1.4f$ for $\zeta_0 = %1.1f$'], eigs(1,i), eigs(2,i), z_0(i))
end
w = warning ('on','all');


