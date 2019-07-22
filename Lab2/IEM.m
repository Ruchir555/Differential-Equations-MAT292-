%Objective: Write your own ODE solver (using the Heun/Improved Euler Method).
%Details: This m-file should be a function which accepts as variables (t0,tN,y0,h), where t0 and tN are the start and end points of the interval 
%on which to solve the ODE, y0 is the initial condition of the ODE, and h is the stepsize. 
%You may also want to pass the function into the ODE the way ode45 does (check lab 2).
%Note: you will need to use a loop to do this exercise. You will also need to recall the Heun/Improved Euler algorithm learned in lectures. 

function [t,y] = IEM(f, t0, tN, y0, h)
   
N = (tN - t0)/ h;
t = linspace(t0, tN, (N+1));
y = zeros(1, length(t));
y(1) = y0;       % Matlab indexes the first element from '1'

for i = 1 : N
    y_1 = f(t(i),y(i));
    y(i+1) = y(i) + h * y_1;
    
    y_2 = f(t(i+1), y(i+1));
    y(i+1) = y(i) + (h/2) * (y_1 + y_2);
end

end