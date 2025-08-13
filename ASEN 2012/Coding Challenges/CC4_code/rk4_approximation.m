function [xout,yout,N] = rk4_approximation(g,del,x0,y0,xf)

N = (xf-x0)/del;

xout = linspace(x0,xf,N);

yout = zeros(N,1);

yout(1)=y0;


for i = 1:(N-1)
    ydot1 = g(xout(i),yout(i));
    ydot2 = g((xout(i) + del/2), (yout(i) + (del/2) * ydot1));
    ydot3 = g((xout(i) + del/2), (yout(i) + (del/2) * ydot2));
    ydot4 = g((xout(i) + del), (yout(i) + del * ydot3));

    yout(i+1) = yout(i) + (del/6)*(ydot1 + 2*ydot2 + 2*ydot3 + ydot4);


end

    
end