%% ES155 P7
%% Problem 3

%% 3.a
P = tf(1, [1 10 3 10])
S = 1000 * tf([1 1], [1 10])

L = P*S

figure(1); clf;
subplot(1,2,1)
bode(L)
subplot(1,2,2)
nyquist(L)

saveas(gcf, "ES155P7_3a.jpg")

pole(L)
[GainMargin, PhaseMargin, Wcg, Wcp] = margin(L)

%% 3.b
P = tf(100, [100, 101, 1])
S = tf([1 10], 1)

L = P*S

figure(2); clf;
subplot(1,2,1)
bode(L)
subplot(1,2,2)
nyquist(L)

saveas(gcf, "ES155P7_3b.jpg")

pole(L)
[GainMargin, PhaseMargin, Wcg, Wcp] = margin(L)

%% Problem 4

% Constants
a = 0.2
b = 25
c = 50
T = 200
m = 1000

P = tf(T*b*a, [m, a*m + c, a*c])

%% 4.a

table = zeros(4,8);

figure(3); clf;
figure(4); clf;

Kp = [0.5, 0.05, 0.05, 0.005]
Ki = [0.1, 1, 0.001, 0.001]

for i = 1:4
   figure(3)
   kp = Kp(i)
   ki = Ki(i)
   
   table(i,1) = kp;
   table(i,2) = ki;
   
   C = tf([kp, ki], 1)
   G = feedback(P*C, 1)
   
   L = P*C
   sysL = ss(L)
   
   sys = ss(G);
   
   table(i,3) = isstable(sys);
   
   [GainMargin, PhaseMargin, Wcg, Wcp] = margin(sysL)
   
   table(i,4) = GainMargin;
   table(i,5) = PhaseMargin;
      
   [y, t] = step(sys);
   if y(length(t)) == y(length(t-2))
       error_SS = 1 - y(length(t))
   else
       error("Steady State not reached")
   end
   
   table(i,6) = error_SS;
   
   S = stepinfo(sys)

   table(i,7) = S.RiseTime;
   table(i,8) = S.Overshoot;
   
   subplot(2,4,i)
   pzmap(sys)
   
   subplot(2,4,i+4)
   step(sys, S.RiseTime*10)
      
   figure(4)
   subplot(2,4,i)
   bode(sys)
   subplot(2,4,i+4)
   nyquist(sys)
   
   % plot the unit circle on the nyquist plot
   theta = 0:0.1:2*pi;
   x = cos(theta);
   y = sin(theta);
   hold on;
   plot(x,y,'--');
   hold off;
end

table







