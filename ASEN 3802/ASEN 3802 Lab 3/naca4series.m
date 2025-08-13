function coords = naca4series(digits,c,n)
    %interpreting the digits    
    t = digits(3)/100; % Percent max thickness
    p = digits(2)/10; % 1/10 loc of max camber
    m = digits(1)/100; % Percent max camber

    % Creating panel nodes
    x = linspace(0,c,n)'; %Equispaced nodes

    % Chebyshev nodes
    for k = 1:n
        x(k) = c/2*cos((2*k-1)*pi/(2*n)) + c/2;
    end

    % Thickness Distribution
    yt = t/(0.2)*c.*(0.2969.*sqrt(x/c)-0.1260.*(x/c)-0.3516.*(x/c).^2+0.2843*(x/c).^3-0.1036.*(x/c).^4);

    % Camber line distribution
    yc = zeros(length(x),1);
    discontinuity = find(x>p*c,1);
    for i = 1:discontinuity
        yc(i) = m*x(i)/p^2*(2*p-x(i)/c);
    end
    for i = discontinuity:length(x)
        yc(i) = m*(c-x(i))/(1-p)^2*(1+x(i)/c-2*p);
    end

    % Camber line derivative
    ycp = zeros(length(x),1);
    for i = 1:discontinuity
        ycp(i) = m*2*p/p^2 - -2*m*x(i)/(p^2*c);
    end
    for i = discontinuity:length(x)
        ycp(i) = m*(-x(i))/(1-p)^2*(1+x(i)/c-2*p) + m*(c-x(i))/(1-p)^2*(1/c);
    end

    % xi
    xi = atan(ycp);

    %Calculating  the upper and lower surfaces
    coords.xU = x - yt.*sin(xi);
    %coords.xU(end) = [];
    coords.xL = x + yt.*sin(xi);
    %coords.xL(1) = [];
    coords.yU = yc + yt.*cos(xi);
    %coords.yU(end) = [];
    coords.yL = yc - yt.*cos(xi);
    %coords.yL(1) = [];
    coords.yc = yc;
    coords.xc = x;
    coords.xb = cat(1,coords.xL,flip(coords.xU));
    coords.yb = cat(1,coords.yL,flip(coords.yU));
    
end