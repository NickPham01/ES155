%% ES155 HW6 Problem 2

m = 1000;
c = 50;
b = 25;
a = 0.2;
T = 200;

%% Part b
kp = [0.01 0.1]

figure(1); clf;
for i = 1:length(kp)
    
    H_eu = tf(kp(i))
    H_uT = tf(a*T, [1, a])
    H_Tv = tf(b, [m, c])
    H_ev = H_Tv * H_uT * H_eu
    H_fdbk = tf(1)
    H = feedback(H_ev, H_fdbk)
    
    sys = ss(H)
    subplot(length(kp),2,2*i-1)
    step(sys)
    
    subplot(length(kp),2, 2*i)
    bode(sys)
    
end

%% Part d

kp = 0.5
ki = 0.1

H_eu = tf([kp, ki], [1, 0])
H_uT = tf(a*T, [1, a])
H_Tv = tf(b, [m, c])
H_fdbk = tf(1, 1)
H = feedback(H_Tv * H_uT * H_eu, H_fdbk)

sys = ss(H)
figure(2); clf;

subplot(1, 2, 1)
step(sys)

subplot(1, 2, 2)
bode(sys)

%% Try using separate transfer functions:

H_eu = tf(kp, 1)
H_uT = tf(a*T, [1, a])
H_Tv = tf(b, [m, c])
H_fdbk = tf(1, 1)
H = feedback(H_Tv * H_uT * H_eu, H_fdbk)

%[A, B, C, D] = tf2ss(H)
sys = ss(H)
figure(3); clf;
subplot(1,2,1)
step(sys)
subplot(1, 2, 2)
bode(sys)


