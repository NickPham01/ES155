%% Problem 3

clear; clc;

%% recreate plot B.1
k = 0.001
a = 10
b = 0.1

B1 = tf(k, [1, 0]) * tf([1, a], [1, b])

figure(1); clf;
subplot(1,3,1)
h1 = bodeplot(B1)
setoptions(h1, 'Grid', 'on')

subplot(1,3,2)
nyquist(B1)

sys1 = feedback(B1, 1)
subplot(1,3,3)
step(sys1)


%% recreate plot B.2

k = 200
w0 = 2
z0 = 0.05
a = 20

B2 = k * tf(1, [1, 2*z0*w0, w0^2]) * tf(1, [1, a])

figure(2);clf;
subplot(1,3,1)
h2 = bodeplot(B2)
setoptions(h2, 'Grid', 'on')

subplot(1,3,2)
nyquist(B2)

subplot(1,3,3)
sys2 = feedback(B2, 1)
step(sys2)

%% recreate plot B.3

a = 20
b = 0.05
c = 1

k = b*c/a * 10


B3 = k * tf([1, a],1) * tf(1, [1, b]) * tf(1, [1, c])

figure(3);clf;
subplot(1,3,1)
h3 = bodeplot(B3)
setoptions(h3, 'Grid', 'on')

subplot(1,3,2)
nyquist(B3)

sys3 = feedback(B3, 1)
subplot(1,3,3)
step(sys3)

