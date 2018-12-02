function showBode(figNum, sys, SSerr, PM, BWerr, BW, titleChars)
figure(figNum); clf;

% get the bode plot values
[mag, phase, wout] = bode(sys);
mag = reshape(mag, [length(mag), 1]);
phase = reshape(phase, [length(phase), 1]);


% make the bode plot
subplot(2,1,1)
loglog(wout, mag)
title({titleChars, '', 'Magnitude'})
ylabel('Magnitude (dB)')

subplot(2,1,2)
semilogx(wout, phase)
title("Phase")
xlabel('\omega (rad/s)')
ylabel('Phase (deg)')

% given the desired tracing error, compute the necessary gain in the BW
% T = L(s)/(1 + L(s)) <= BWerr, so L(s) > 1/BWerr in the BW
BWgain = 1/BWerr

subplot(2,1,1)
Xlim = xlim;
Ylim = ylim;
x = Xlim(1);
y = Ylim(1);
w = BW - x;
h = BWgain - y;
rectangle('Position', [x, y, w, h])

% show the desired phase margin
minPhaseAtBW = -180 + PM;
[Gm, Pm, Wcg, Wcp] = margin(sys)

subplot(2,1,2)
%{
Xlim = xlim;
Ylim = ylim;
x = Xlim(1);
y = Ylim(1);
w = Wcg - x;
h = minPhaseAtBW - y;
rectangle('Position', [x,y,w,h])
%}
figure(figNum+1); clf;
margin(sys)

end

