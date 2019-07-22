function [tc, yc] = adaptive_euler1(f, t0 , tN, y0, h)
   
%  |{ [ f(t[n+1],y[n]) - f(t[n],y[n]) ] + [ f(t[n],y[n+1]) - f(t[n],y[n]) ] } / h|
t = t0;
tc = [t0];
y = y0;
yc = [y0];
Cprev = ((f(t+h,y(1)) - f(t(1),y(1)) ) + ( f(t(1),y(1)+h*f(t(1),y(1))) - f(t(1),y(1)))) / h;
%disp(Cprev);

while(t < tN)
    hvar = h;  
    k1 = f(t,y);
    tnext = t + hvar;
    ynext = y + hvar * k1;
    C = (f(tnext, y) - f(t, y) + f(t, ynext) - f(t, y)) / hvar;

    while (C > 1.1*Cprev & hvar > 0.0001)
        hvar = hvar/2;
        tnext = t + hvar;
        ynext = y + hvar * k1;
        C = (f(tnext, y) - f(t, y) + f(t, ynext) - f(t, y)) / hvar;
        %disp(C);
    end
         
    t = tnext;
    tc = [tc, t];
    y = ynext;
    yc = [yc, y];
    Cprev = C;
end
end