function [data] = nc_read(nc_filename)

    nc = netcdf.open(nc_filename, 'NC_NOWRITE');
    
    
    thk_id = netcdf.inqVarID(nc, 'thk');
    thk_t_id = netcdf.inqVarID(nc,'thk_t');
    %stagthk_id = netcdf.inqVarID(nc,'stagthk');
    uvel_id = netcdf.inqVarID(nc, 'uvelhom');
    vvel_id = netcdf.inqVarID(nc, 'vvelhom');
    time_id = netcdf.inqVarID(nc, 'time');
    bmlt_id = netcdf.inqVarID(nc, 'bmlt');
    
    x0_id = netcdf.inqVarID(nc, 'x0');
    y0_id = netcdf.inqVarID(nc, 'y0');
    x1_id = netcdf.inqVarID(nc, 'x1');
    y1_id = netcdf.inqVarID(nc, 'y1');
    
    %time = netcdf.getVar(nc,time_id);
    fstagx = @(x) x(:);
    fstagy = @(y) y(1:(end-3));
    fgridx = @(x) x(2:(end-1));
    fgridy = @(y) y(2:(end-4));
    %f1 = @(M) M(2:(end-1),2:(end-4),:);
    fgrid = @(M) M(2:(end-1),2:(end-4),:,:);
    fstag = @(M) M(:,1:(end-3),:);
    fstag_layer = @(M) M(:,1:(end-3),:,:);
    data.xstag = fstagx(netcdf.getVar(nc, x0_id));
    data.xgrid = fgridx(netcdf.getVar(nc, x1_id));
    data.ystag = fstagy(netcdf.getVar(nc, y0_id));
    data.ygrid = fgridy(netcdf.getVar(nc, y1_id));
    data.thk = fgrid(netcdf.getVar(nc, thk_id));
    data.thk_t = fgrid(netcdf.getVar(nc,thk_t_id));
    %data.stagthk= fstag(netcdf.getVar(nc,stagthk_id));
    data.uvel = fstag_layer(netcdf.getVar(nc, uvel_id));
    data.vvel = fstag_layer(netcdf.getVar(nc, vvel_id));
    data.bmlt = fgrid(netcdf.getVar(nc,bmlt_id));
    
    [m,n,k] = size(data.thk);
    
    
    data.uvelmean = squeeze(mean(data.uvel,3));
    data.vvelmean = squeeze(mean(data.vvel,3));
    %data.uflx = data.uvelmean.*data.stagthk;
    %data.vflx = data.vvelmean.*data.stagthk;
    
    dx = data.xstag(2)-data.xstag(1);
    dy = data.ystag(2)-data.ystag(1);
    
    ux = (data.uvelmean(2:end,    1:(end-1),:) + ...
          data.uvelmean(2:end,    2:end,:)- ...
          data.uvelmean(1:(end-1),1:(end-1),:)- ...
          data.uvelmean(1:(end-1),2:end,:))/(2*dx);
    vy = (data.vvelmean(1:(end-1),2:end,:)+ ...
          data.vvelmean(2:end,    2:end,:)- ...
          data.vvelmean(1:(end-1),1:(end-1),:)- ...
          data.vvelmean(2:end,    1:(end-1),:))/(2*dy);
      
  
    data.vel_div = (ux+vy) .* data.thk;
    
    %data.flux_div = ...
    %     (data.uflx(2:end,    1:(end-1),:) + ...
    %      data.uflx(2:end,    2:end,:)- ...
    %      data.uflx(1:(end-1),1:(end-1),:)- ...
    %      data.uflx(1:(end-1),2:end,:))/(2*dx) + ...
    %     (data.vflx(1:(end-1),2:end,:)+ ...
    %      data.vflx(2:end,    2:end,:)- ...
    %      data.vflx(1:(end-1),1:(end-1),:)- ...
    %      data.vflx(2:end,    1:(end-1),:))/(2*dy);

    %data.y_adv = zeros(m,n,k);
    %data.y_adv(:,2:(end-1),:) = ...
    %    data.vvelmean(round(m/2),1,1) * ...
    %    (data.thk(:,3:end,:) - data.thk(:,1:(end-2),:))/dy;
        
    [X1,Y1] = meshgrid(data.xgrid,data.ygrid);
    data.Xgrid = X1';
    data.Ygrid = Y1';
    

end