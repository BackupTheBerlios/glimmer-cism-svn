function [ x,u,h,w ] = steady_ice_1(A,rhoi,rhoo,g,x0,x1,n, u0, h0, a_net_per_year)
%steady_ice_1 compute 1d ice profile with imposed upstream flux


if (a_net_per_year ~= 0.0) 
    wfun = @(x) a_net_per_year*(x-x0)+u0*h0;
    ufun = @(x) u0*(1 + (A/a_net_per_year)*((1-rhoi/rhoo)*rhoi*g/4)^3*((wfun(x)/u0).^4 - (h0)^4)).^(1/4);

    x = x0:(x1-x0)/n:x1;
    u = ufun(x);
    h = wfun(x)./u;
    w = wfun(x);

else
    wfun = @(x) u0*h0;
    ufun = @(x) u0*(1 + (4*A)*((1-rhoi/rhoo)*rhoi*g*h0/4)^3*(x-x0)/u0).^(1/4);

    x = x0:(x1-x0)/n:x1;
    u = ufun(x);
    h = wfun(x)./u;
    w = wfun(x);
end 

end

