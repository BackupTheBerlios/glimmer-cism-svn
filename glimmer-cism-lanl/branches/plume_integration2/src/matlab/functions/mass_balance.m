function [melt,in_flux,out_flux,total_acab,unsteady] = mass_balance(dice)

[m,n,k] = size(dice.bmlt);

tslice = k;
j_front = n;
j_inflow = 1;

dx = dice.xstag(2)-dice.xstag(1);
dy = dice.ystag(2)-dice.ystag(1);
dt = dice.time(end)-dice.time(end-1);

width = dice.xstag(end)-dice.xstag(1);
len = dice.ystag(end)-dice.ystag(1);

gradx = zeros(m,n);
grady = zeros(m,n);
gradx(2:end-1,:) = (dice.thk(3:end,:,tslice)-dice.thk(1:end-2,:,tslice))/(2*dx);
gradx(1,:)       = (dice.thk(2    ,:,tslice)-dice.thk(1,:,tslice))/dx;
gradx(end,:)     = (dice.thk(end  ,:,tslice)-dice.thk(end-1,:,tslice))/dx;
grady(:,2:end-1) =  (dice.thk(:,3:end,tslice)-dice.thk(:,1:end-2,tslice))/(2*dy);
grady(:,1)       = (dice.thk(:,2    ,tslice)-dice.thk(:,1,tslice))/dy;
grady(:,end)     = (dice.thk(:,end  ,tslice)-dice.thk(:,end-1,tslice))/dy;
grad_norm = sqrt(1+gradx.^2+grady.^2);

acab = -1.0;
total_acab = acab*width*len;

melt = sum(sum(dice.bmlt(:,:,tslice).*grad_norm))*dx*dy;
melt = sum(sum(dice.bmlt(:,:,tslice)))*dx*dy;

%in_flux = sum(0.25*(dice.vvelmean(1:end-1,j_inflow,tslice)+ ...
%                     dice.vvelmean(2:end  ,j_inflow,tslice)+ ...
%                     dice.vvelmean(1:end-1,j_inflow+1,tslice)+...
%		     dice.vvelmean(2:end  ,j_inflow+1,tslice)) .* ...
%	       dice.thk(:,j_inflow,tslice)) * dx ;

in_flux = sum(0.5*(dice.vvelmean(1:end-1,j_inflow,tslice)+ ...
                   dice.vvelmean(2:end  ,j_inflow,tslice)).* ...
	      dice.thk(:,j_inflow,tslice)) * dx ;


out_flux = sum(0.25*(dice.vvelmean(1:end-1,j_front,tslice)+ ...
                     dice.vvelmean(2:end  ,j_front,tslice)+ ...
                     dice.vvelmean(1:end-1,j_front-1,tslice)+...
		     dice.vvelmean(2:end  ,j_front-1,tslice)) .* ...
	       dice.thk(:,j_front,tslice)) * dx ;




unsteady = sum(sum((dice.thk(:,:,end)-dice.thk(:,:,end-1))/dt * (dx*dy)));
end
