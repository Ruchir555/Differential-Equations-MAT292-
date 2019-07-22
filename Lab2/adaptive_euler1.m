function [tc, yc] = adaptive_euler1(f, t0 , tN, y0, h)
   
%  |{ [ f(t[n+1],y[n]) - f(t[n],y[n]) ] + [ f(t[n],y[n+1]) - f(t[n],y[n]) ] } / h|
t = t0;
tc = [t0];
y = y0;
yc = [y0];
Cbefore = ((f(t+h,y(1)) - f(t(1),y(1)) ) + ( f(t(1),y(1)+h*f(t(1),y(1))) - f(t(1),y(1)))) / h;
%disp(Cprev);

while(t < tN)
    h_var = h;  
    k_1 = f(t,y);
    t_next = t + h_var;
    y_next = y + h_var * k_1;
    C = (f(t_next, y) - f(t, y) + f(t, y_next) - f(t, y)) / h_var;

    while (C > 1.1*Cbefore & h_var > 0.0001)
        h_var = h_var/2;
        t_next = t + h_var;
        y_next = y + h_var * k_1;
        C = (f(t_next, y) - f(t, y) + f(t, y_next) - f(t, y)) / h_var;
        %disp(C);
    end
         
    t = t_next;
    tc = [tc, t];
    y = y_next;
    yc = [yc, y];
    Cbefore = C;
end
end