function [y_sol,thk_sol,vel_sol, thk_paterson] = quasi_exact_shear_sol(Vin,Hin,yin,yout,acab,tauxy0,L,ny)

%define constants
A = 10.0^(-16);
rhoi = 910.0;
rhoo = 1028.0;
g = 9.81;

% constants defined for convenience
k = (1-rhoi/rhoo)*(rhoi*g/4.0);
q = tauxy0/(2*L);

% The equation being solved is:

%  dv/dy = A * (k * H - q*w/H) ^3
%  dw/dy = H         
%  d(vH)/dy = acab

% where w = int_{y0}^y H dy is the area of contact with 
% the side walls from the ice front back to the current position and 
% where A,k,q,acab are given constants.

% We can solve for H immediately in terms of v, yielding

% dv/dy = A* ( k* flux(y)/v - q*w*v/flux(y) ) ^ 3
% dw/dy = flux(y)/v

% where flux(y) is defined by:
flux = @(y) Vin*Hin - acab*(yin - y);

% If we define the unknown by u where
% u = [v; w]
% then the ode is written 
% u' = odefun(y,u) where:
odefun = @(y,u) [ A*(k*flux(y)/u(1) - q*u(2)*u(1)/flux(y))^3.0;
                  flux(y)/u(1) ];

% The boundary conditions we use are
% w(icefront) = 0
% v(inflow) = vin

% This can be written in terms of u as:
bcfun = @(u0, u1) [u0(2); u1(1)-Vin];

% we need to define an initial guess at a solution
dy = (yin-yout) / (ny-1);  % step size of numerical solution 
solinit.x = yout:dy:yin;   % discrete values of independent
                           % variable y (needs to be called 'x')
solinit.y = zeros(2,ny);  
solinit.y(1,:) = Vin;      % Guess a uniform velocity equal to the
                           % inflow velocity
solinit.y(2,:) = 0:(Hin*yin*0.5/(ny-1)):(Hin*yin*0.5); % Guess for w
                                                     % corresponding to
                                                     % constant v
                                                     % and H
                                                     % increasing linearly

% ask matlab to solve the bvp
sol = bvp4c(odefun, bcfun,solinit);

y_sol = sol.x;
vel_sol = sol.y(1,:);  % extract the first component of the sol
                     % Note: we don't care about w = u(2,:)
thk_sol = flux(y_sol) ./ vel_sol;  % calculate thickness solution

paterson_fun = @(y) Hin + 2*tauxy0/(rhoi*g*L*(1-rhoi/rhoo))*(y-yin);

thk_paterson = paterson_fun( y_sol );
end
