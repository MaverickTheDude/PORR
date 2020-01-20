function [y, z] = diff24(x,t)
%DIFF24 is a second derivative 4th order of accuracy finite differences
%formula for a vector x. Constant step within each row of x is assumed

% Number of points which will be used to interpolate four missing points after deriviation.
points = 3;
[m, n] = size(x);
y = zeros(m,n);
for j = 1:m
    dt = t(2) - t(1);
    for i = 3:n-2
        y(j,i) = (-x(j,i-2) + 16*x(j,i-1) - 30*x(j,i) + 16*x(j,i+1) - x(j,i+2) ) / (12*dt^2);
    end
    
    % Interpolation of two initial and two last points
    order = points-1;
    sampleBeg = y(j,3:points+2);
    sampleEnd = y(j,end-2-points+1:end-2);
    domainBeg = t(3:points+2);
    domainEnd = t(end-2-points+1:end-2);
    p = polyfit(domainBeg,sampleBeg,order);
    y1 = polyval(p,t(1:2));
    p = polyfit(domainEnd,sampleEnd,order);
    y2 = polyval(p,t(end-1:end));
    y(j,1:2) = y1; y(j,end-1:end) = y2;
end

% Second output is the first derivative 2nd order accurate
if nargout == 2
z = zeros(m,n);
for j = 1:m
    for i = 3:n-2
        z(j,i) = (x(j,i-2) - 8*x(j,i-1) + 8*x(j,i+1) - x(j,i+2) ) / (12*dt);
    end
    
    % Interpolation of two initial and two last points
    order = points-1;
    sampleBeg = z(j,3:points+2);
    sampleEnd = z(j,end-2-points+1:end-2);
    domainBeg = t(3:points+2);
    domainEnd = t(end-2-points+1:end-2);
    p = polyfit(domainBeg,sampleBeg,order);
    y1 = polyval(p,t(1:2));
    p = polyfit(domainEnd,sampleEnd,order);
    y2 = polyval(p,t(end-1:end));
    z(j,1:2) = y1; z(j,end-1:end) = y2;
end 
end