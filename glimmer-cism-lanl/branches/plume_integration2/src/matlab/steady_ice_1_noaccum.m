
function [ x,u,h,w ] = steady_ice_1_noaccum(A,rhoo,rhoi,g,x0,x1,n, u0, h0)
%steady_ice_1 compute 1d ice profile with imposed upstream flux
%and no accumulation

wfun = @(x) u0*h0;
ufun = @(x) u0*(1 + (4*A)*((1-rhoi/rhoo)*rhoi*g*h0/4)^3*(x-x0)/u0).^(1/4);

x = x0:(x1-x0)/n:x1;
u = ufun(x);
h = wfun(x)./u;
w = wfun(x);

end


