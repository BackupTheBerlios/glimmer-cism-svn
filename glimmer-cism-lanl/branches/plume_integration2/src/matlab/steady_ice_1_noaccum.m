
function [ x,u,h,w ] = steady_ice_1_noaccum( x0,x1,n, u0, h0)
%steady_ice_1 compute 1d ice profile with imposed upstream flux
%and no accumulation

rho = 930;
g = 9.8;
A = 1*10^(-20);

wfun = @(x) u0*h0;
ufun = @(x) u0*(1 + (4*A/u0)*(rho*g*h0/4)^3*(x-x0)).^(1/4);

x = x0:(x1-x0)/n:x1;
u = ufun(x);
h = wfun(x)./u;
w = wfun(x);

end


