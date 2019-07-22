%Objective: Create an Adaptive Euler method, with an adaptive step size h.
%Details: Create an m-file which accepts the variables (t0,tN,y0,h), as in exercise 1, where h is an initial step size. 
%You may also want to pass the function into the ODE the way ode45 does.
%Create an implementation of Euler's method by modifying your solution to exercise 1. Change it to include the following:
%(a) On each timestep, make two estimates of the value of the solution at the end of the timestep:
%Y from one Euler step of size h and Z from two successive Euler steps of size h/2. The difference in these two values is an estimate for the error.
%(b) Let tol=1e-8 and D=Z-Y. If abs(D)<tol, declare the step to be successful and set the new solution value to be Z+D. This value has local error O(h^3). 
%If abs(D)>=tol, reject this step and repeat it with a new step size, from (c).
%(c) Update the step size as h = 0.9*h*min(max(tol/abs(D),0.3),2).
%Comment on what the formula for updating the step size is attempting to achieve.

function [tc,yc] = adaptive_euler(f, t0, tN, y0, h)
   
N = (tN - t0)/h;
tc = linspace(t0, tN, (N+1));
yc = zeros(1, length(tc));
y=y0;
t=t0;
tc = [t0];
yc = [y0];       % Matlab indexes the first element from '1'
tol = 10^(-8);

   while (t<tN) 
       h_var = h;
       k1 =f(t,y);
       t_next = t + h_var;
       y_next = y + h_var * k1;
       y_1 = y + (h_var/2) * k1;
       y_2 = y + (h_var/2) * y_1;
       error = y_next - y_2;
        
    while (abs(error))>=tol
        h_var = 0.9*h*min(max(tol/abs(error),0.3),2);
        t_next = t + h_var;
        y_next = y + h_var * k1;
        y_1 = y + (h_var/2) * k1;
        y_2 = y + (h_var/2) * y_1;
        error = y_next - y_2;
    end
    
        t = t_next;
        tc = [tc, t];
        yc = [yc, y];
        y = error + y_2;  %Z+D
   end         
end
    
















