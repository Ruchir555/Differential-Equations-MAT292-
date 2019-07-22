function [t,y] = DE2_tulluruc(p, q, g, t0, tN, y0, y1, h)

    %Solves a second order differential equation of the form:
    %                y'' + p(t)y' + q(t)y = g(t)

    N = floor((tN - t0)/h); % number of steps

    %Initialize arrays
    t = zeros(1, N+1);
    y = zeros(1, N+1);
    y_prime = zeros(1, N+1);

    %Initial conditions
    t(1) = t0;
    t(2) = t0 + h;

    y(1) = y0;
    y(2) = y0 + h*y1;

    y_prime(2) = (y(2) - y(1))/h;

    %Loop
    for i = 2:N
        y(i+1) = (h^2)*(-p(t(i))*y_prime(i) - q(t(i))*y(i) + g(t(i))) + 2*y(i) - y(i-1);    %Rearranged finite differences method
        y_prime(i+1) = (y(i+1) - y(i))/h;   %Y_prime calculation based on y(i+1) calculated in this step
        t(i+1) = t(i) +h;   %Update timestep
    end
end

    
    
    
    
    
    
    
    
    
    
    
