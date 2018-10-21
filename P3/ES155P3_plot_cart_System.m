function ES155P3_plot_cart_System(ic, fignum)
tspan = [0 15]

 % compute output
[t, x] = ode45(@(t,x) ES155P3_cart_inv_pend(t, x), tspan, ic);

size(x)
figure(fignum);clf;

subplot(3,1,1)
plot(t, x(:,2))
ylabel('$\theta$', "Interpreter", "latex")
title('$\theta$', "Interpreter", "latex")
xlabel("t", 'interpreter', 'latex')

subplot(3,1,2)
plot(t, x(:,3))
ylabel('$\dot{p}$', "Interpreter", "latex")
title('$\dot{p}$', "Interpreter", "latex")
xlabel("t", 'interpreter', 'latex')


subplot(3,1,3)
plot(t, x(:,4))
ylabel('$\dot{\theta}$', "Interpreter", "latex")
title('$\dot{\theta}$', "Interpreter", "latex")
xlabel("t", 'interpreter', 'latex')

annotation('textbox', [0 0.9 1 0.1], ...
    'String', 'Cart Inverted Pendulum Sysetm', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center')

end