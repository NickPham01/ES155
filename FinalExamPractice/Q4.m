P = tf(50, [1, 0.1, 0])
C = tf(1)

L = P*C

bode(L)

[Gm,Pm,Wcg,Wcp] = margin(L)
