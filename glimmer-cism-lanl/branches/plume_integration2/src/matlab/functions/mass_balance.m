function [melt,applied_melt,in_flux,out_flux,total_acab,dvdt] = mass_balance(dice,dplume,tslices)

[m,n,k] = size(dice.bmlt);

if (tslices < 0)
  tslices = [k];
end

j_front = n;
j_inflow = 1;

dx = dice.xstag(2)-dice.xstag(1);
dy = dice.ystag(2)-dice.ystag(1);
if (length(dice.time) < 2) 
  error('need at least two times to determine dt'); %dt = 0.1;
else
  dt = dice.time(end)-dice.time(end-1);
end

width = dice.xstag(end)-dice.xstag(1);
len = dice.ystag(end)-dice.ystag(1);

gradx = zeros(m,n);
grady = zeros(m,n);

melt = zeros(length(tslices),1);
applied_melt = zeros(length(tslices),1);
in_flux = zeros(length(tslices),1);
out_flux = zeros(length(tslices),1);
dvdt = zeros(length(tslices),1);


for i=1:length(tslices)
  tslice = tslices(i);

gradx(2:end-1,:) = (dice.thk(3:end,:,tslice)-dice.thk(1:end-2,:,tslice))/(2*dx);
gradx(1,:)       = (dice.thk(2    ,:,tslice)-dice.thk(1,:,tslice))/dx;
gradx(end,:)     = (dice.thk(end  ,:,tslice)-dice.thk(end-1,:,tslice))/dx;
grady(:,2:end-1) =  (dice.thk(:,3:end,tslice)-dice.thk(:,1:end-2,tslice))/(2*dy);
grady(:,1)       = (dice.thk(:,2    ,tslice)-dice.thk(:,1,tslice))/dy;
grady(:,end)     = (dice.thk(:,end  ,tslice)-dice.thk(:,end-1,tslice))/dy;
grad_norm = sqrt(1+gradx.^2+grady.^2);

acab = -1.0;
total_acab = acab*width*len;

applied_melt(i) = sum(sum(dice.bmlt(:,:,tslice)))*dx*dy;

if(dplume.has_bmelt_avg)
  melt(i) = sum(sum(dplume.bmelt_avg(:,:,tslice)))*dx*dy;
    if (mean(mean(dplume.bmelt_avg(:,:,tslice))) > 1000.0)
      %must be a badly scaled case so fix it
	 melt(i) = melt/(3600.0*24.0*365.25);
    end
else
  melt(i) = sum(sum(dice.bmlt(:,:,tslice)))*dx*dy;
end

in_flux(i) = sum(0.5*(dice.vvelmean(1:end-1,j_inflow,tslice)+ ...
                   dice.vvelmean(2:end  ,j_inflow,tslice)).* ...
	      dice.thk(:,j_inflow,tslice)) * dx ;

out_flux(i) = sum(0.25*(dice.vvelmean(1:end-1,j_front,tslice)+ ...
                     dice.vvelmean(2:end  ,j_front,tslice)+ ...
                     dice.vvelmean(1:end-1,j_front-1,tslice)+...
		     dice.vvelmean(2:end  ,j_front-1,tslice)) .* ...
	       dice.thk(:,j_front,tslice)) * dx ;


dvdt(i) = sum(sum((dice.thk(:,:,tslice)-dice.thk(:,:,tslice-1))/dt * (dx*dy)));

end

end
