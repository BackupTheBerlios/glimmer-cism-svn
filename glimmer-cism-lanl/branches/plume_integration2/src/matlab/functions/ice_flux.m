function [total_flux] = ice_flux(dice)


  vvel = dice.vvelmean(:,:,end);
grid_vvel = 0.25*(vvel(1:end-1,end)+vvel(2:end,end)+vvel(1:end-1,end-1)+vvel(2:end,end));
total_flux = sum(grid_vvel .* dice.thk(:,end,end)) * (dice.xstag(2)-dice.xstag(1))/1e9;

end
