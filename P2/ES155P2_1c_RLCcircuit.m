function dydt = ES155P2_1c_RLCcircuit(t, y, omega_0, gains)
    % Define Constants
    R = 1;
    L = 0.1;
    C = 0.2;

    % calculate controller input based on time (unit step at t=0)
    u = 1;
    
    dydt = zeros(2,1);
    dydt(1) = (y(2) - y(1)/R)/C;
    dydt(2) = (u - y(1))/L;
end