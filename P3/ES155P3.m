%% ES155 P3


%% 1.c

A_100_20 = [0.1 -0.5; 20*20/200^2 (20/200 - 0.1)]
eig(A_100_20)

A_0_0 = [0.1 -0.5; 20/100^2 -0.1]
eig(A_0_0)

%% 1.d

omega = -0.61
A_bar = [0.1+omega -0.5; 0.01 0]
eig(A_bar)


%% 2.c

M = 10
m = 80
I = 100
l = 1
g = 9.8
c = 0.1
gamma = 0.01

denom = (M + m)*(I + m*l^2) - m^2*l^2

A_0 = [ 0, 0, 1, 0;
        0, 0, 0, 1;
        0, m^2*l^2*g/denom, -c*(I + m*l^2)/denom, -gamma*l*m/denom;
        0, (M + m)*m*g*l/denom, -c*l*m/denom, -gamma*(M + m)/denom]

latex(vpa(sym(A_0), 3))

    
A_pi = A_0 .* [1 1 1 1;
             1 1 1 1;
             1 1 1 -1;
             1 -1 -1 1]
         
latex(vpa(sym(A_pi), 3))

eig(A_0)
eig(A_pi)

%% 2.d

B = [0; 0; (I + m*l^2)/denom; l*m/denom]
K = [-15.3 1730 -50 443]

B*K

A_bar = A_0 + B*-K

eig(A_bar)

latex(vpa(sym(A_bar),3))

%% 2.e

ic = [0;1;0;0]
ES155P3_plot_cart_System(ic,1)
saveas(gca, "UnstableIC_0_1_0_0.png")


ic = [0; 0.5; 0; 0]
ES155P3_plot_cart_System(ic,2)
saveas(gca, "StableIC_0_05_0_0.png")




















