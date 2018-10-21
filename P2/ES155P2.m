%% Problem 1

%% Part b

maxrangex1 = 5
maxrangex2 = 5
stepsizex1 = 1
stepsizex2 = 1

R = 1;
L = 0.1;
C = 0.2;


figure(1); clf;
[x1, x2] = meshgrid(-maxrangex1:stepsizex1:maxrangex1, -maxrangex2:stepsizex1:maxrangex2);
x1dot = (x2 - x1/R)./C;
x2dot = -x1./L;
quiver(x1,x2,x1dot, x2dot);
xlabel("x1")
ylabel("x2")

saveas(gcf, "ES155P2_1b_phaseportrait.png")

%% Part c with ode45

tspan = [0 5]
ic = [0; 0]         % given in problem

 % compute output
[t, y] = ode45(@(t,y) ES155P2_1c_RLCcircuit(t, y), tspan, ic);

size(y)
figure(2);clf;
plot(t, y)
xlabel("t", 'interpreter', 'latex')
ylabel("value")
legend("Voltage", "Current")
title("$V$ and $I$ Response of RLC circuit", 'interpreter', 'latex')

saveas(gcf, "ES155P2_1c_VIresponse.png")


%% Problem 3
%% Part b

A = [1 1; 1 2]

rref(A)
rank(A)


%% Part d
B = [1 1 0; 1 1 1]

rref(B)
rank(B)

%% Part f

C = B'

rref(C)
rank(C)