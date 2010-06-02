function [ x,u,h,w ] = steady_ice_1( x0,x1,n, u0, h0, a_net_per_year)
%steady_ice_1 compute 1d ice profile with imposed upstream flux


a_net = a_net_per_year / (3600*24*365.25);

rho = 930;
g = 9.8;
A = 5*10^(-18);

wfun = @(x) a_net*(x-x0)+u0*h0;
ufun = @(x) (u0^4 + (A/a_net)*(rho*g/4)^3*(wfun(x).^4 - (u0*h0)^4)).^(1/4);

x = x0:(x1-x0)/n:x1;
u = ufun(x);
h = wfun(x)./u;
w = wfun(x);

end

