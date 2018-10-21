function dxdt = ES155P3_cart_inv_pend(t, x)

    M = 10;
    m = 80;
    l = 1;
    I = 100;
    c = 0.1;
    gamma = 0.01;
    g = 9.8;

    I_t = I + m*l^2;
    M_t = M + m;
    
    dxdt = zeros(4,1);
    
    p = x(1);
    theta = x(2);
    p_dot = x(3);
    theta_dot = x(4);
    
    sine = sin(theta);
    cosine = cos(theta);
    
    K = [-15.3, 1730, -50, 443];
    
    u = -K*x;
    
    dxdt(1) = p_dot;
    dxdt(2) = theta_dot;
    %{
    dxdt(3) = (u - m*l*sine*theta_dot^2 - c*p_dot ...
               - m*l*cosine*gamma/I_t*theta_dot  ...
               + m^2*g*l^2/I_t*cosine*sine)/ ...
               (M_t - m^2*l^2/I_t*cosine^2)
               
    dxdt(4) = (m*l*cosine*u - m*l*cosine*c*p_dot ... 
               - m^2*l^2*cosine*sine*theta_dot^2 ...
               - M_t*(gamma*theta_dot - m*g*l*sine)) / ...
               (M_t * I_t - l^2 * cosine^2)
               %}
           
    dxdt(3) = (-m*l*sine*theta_dot^2 + m*g*(m*l^2/I_t)*sine*cosine - c*p_dot ...
        -gamma/I_t*m*l*cosine*theta_dot + u) / ...
        (M_t - m*(m*l^2/I_t)*cosine^2);
    dxdt(4) = (-m*l^2*sine*cosine*theta_dot^2 + M_t*g*l*sine - c*l*cosine*p_dot ...
        - gamma*M_t/m*theta_dot + l*cosine*u) / ...
        (I_t*M_t/m - m*(l*cosine)^2);
    
    

end