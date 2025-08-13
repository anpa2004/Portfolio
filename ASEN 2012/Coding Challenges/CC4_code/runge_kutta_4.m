function [xout,yout] = runge_kutta_4(g,del,x0,y0,xf)

N = (x0-xf)/del;
xout = linspace(x0,xf,N);

yout(1)=y0;

for i = 1:N
    ydot1 = g(x0,y0);
    ydot2 = g(x0+del/2,y0+(del/2)*ydot1);
    ydot3 = g(x0+del/2,y0+(del/2)*ydot2);
    ydot4 = g(x0+del,y0+del*ydot3);

    yout(i+1) = del/6*(ydot1+2*ydot2+2*ydot3+ydot4);


end

    
end