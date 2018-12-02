%% Part 1
%% 1.a

% Given Constants
g   = 9.81;         % gaccel[m/s^2]
mp  = 0.230;        % massofpendulum[kg]
l   = 0.6413;       % lengthofpendulum[m]
r   = l/2;          % radiustoCOMofpendulum[m]
J   = (1/3)*mp*l^2; % inertiaofpendulumrotatingabout1end[kg-m^2]
y   = 0.0024;       % pendulumdamping[N-m*s]
mc  = 0.38;         % massofcart[kg]
c   = 0.90;         % cartdamping[N-s/m]

% Derived Constants
Mhat = mp + mc;
Jhat = J + mp*r^2;
mu = mp^2 * r^2 - Jhat^2 * Mhat^2;

A = [0 0 1 0;
     0 0 0 1;
     0, g/mu, (Jhat*c)/mu, -(y*mp*r)/mu;
     0, -(Mhat*mp*r*g)/mu, -(mp*r*c)/mu, (Mhat*y)/mu]
 
B = [0; 0; -Jhat/mu; -(mp*r)/mu]

C = [1 0 0 0; 0 1 0 0]

D = [0; 0]

%% 1.b

sys = ss(A, B, C, D)
eig(A)

opt = stepDataOptions('StepAmplitude', 0.1)
[y, t, x] = step(sys, opt);

figure(1); clf;
plot(t, y)
hline = refline(0, 0.1)
hline.Color = 'm'

legend("Position", "Slope", "Input", 'Location', 'Northwest')
saveas(gca, 'ES155Lab2_1b_step.jpg')

%% 1.c


% the smaller poles allow the pendulum to be pushed around before reaching
% an equilibrium, while the large value poles keep the pendulum at the
% equilibrium point theta = 0

figure(2); clf;

plotCount = 1;
pMultipliers = [1, 2, 5, 10]
for i = 1:length(pMultipliers)
    p = [-1, -2, -3, -4];
    p = p.*pMultipliers(i)
    K = place(A, B, p)

    sys = ss(A- B*K, B, C, 0);
    opt = stepDataOptions('StepAmplitude', 0.1);
    [y, t, x] = step(sys, opt);

    titles = ["Position"; "Angle"];
    ylabels = ["$x$", "$\theta$"]
    for j = 1:2
        subplotIdx = plotCount + j -1
        subplot(length(pMultipliers),2, subplotIdx)
        plot(t, y(:,j))
        hline = refline(0, 0.1);
        hline.Color = 'm';
        title({char(titles(j)), ['\lambda = ', num2str(p(1)), ', ', num2str(p(2)), ', ', num2str(p(3)), ', ', num2str(p(4))]})
        ylabel(char(ylabels(j)), 'Interpreter', 'latex')
    end
    
    plotCount = plotCount + 2;
end


subplot(length(pMultipliers), 2, plotCount - 2)
xlabel('$t$', 'Interpreter', 'latex')

subplot(length(pMultipliers), 2, plotCount - 1)
xlabel('$t$', 'Interpreter', 'latex')



