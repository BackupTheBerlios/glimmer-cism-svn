function [ x,u,h,w ] = steady_ice_1(A,rhoi,rhoo,g,x0,x1,n, u0, h0, a_net_per_year)
%steady_ice_1 compute 1d ice profile with imposed upstream flux


if (a_net_per_year ~= 0.0) 
    wfun = @(x) u0*h0 - a_net_per_year*(x0-x);
    ufun4 = @(x) u0^4 - ...
      (A/a_net_per_year)*((1-rhoi/rhoo)*rhoi*g/4)^3 * ...
      ( wfun(x0)^4 - wfun(x).^4);

    x = x0:(x1-x0)/(n-1):x1;
    u = -1.0* ufun4(x) .^(1/4);
    h = wfun(x)./u;
    w = wfun(x);

else
    wfun = @(x) u0*h0;
    ufun4 = @(x) u0^4 - (4*A)*((1-rhoi/rhoo)*rhoi*g/4)^3*(u0*h0)^3*(x0-x);
    x = x0:(x1-x0)/(n-1):x1;
    u = -1.0* ufun4(x) .^(1/4);
    h = u0*h0./u;
    w = wfun(x);
end 

end

