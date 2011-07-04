function [data] = nc_ice_read(nc_filename, timestride, max_slices)

    nc = netcdf.open(nc_filename, 'NC_NOWRITE');
   
    
    thk_id = netcdf.inqVarID(nc, 'thk');
    thk_t_id = netcdf.inqVarID(nc,'thk_t');
       
    lsurf_id = netcdf.inqVarID(nc,'lsurf');
    uvel_id = netcdf.inqVarID(nc, 'uvelhom');
    vvel_id = netcdf.inqVarID(nc, 'vvelhom');
    %time_id = netcdf.inqVarID(nc, 'time');
    bmlt_id = netcdf.inqVarID(nc, 'bmlt');
    x0_id = netcdf.inqVarID(nc, 'x0');
    y0_id = netcdf.inqVarID(nc, 'y0');
    x1_id = netcdf.inqVarID(nc, 'x1');
    y1_id = netcdf.inqVarID(nc, 'y1');
    level_id =netcdf.inqVarID(nc,'level');
    time_dim_id = netcdf.inqDimID(nc,'time');
    time_var_id = netcdf.inqVarID(nc,'time');
    
    try
        uflx_conv_id = netcdf.inqVarID(nc,'uflx_conv');
        vflx_conv_id = netcdf.inqVarID(nc,'vflx_conv');
        have_flux_conv = true;
    catch 
        have_flux_conv = false;
    end
    
    [tname,tlen] = netcdf.inqDim(nc,time_dim_id);    
    
    n_timeslices = min(max_slices,round(floor(tlen/timestride)));
    
    time = netcdf.getVar(nc,time_var_id,0,n_timeslices,timestride);
    
    x0 = netcdf.getVar(nc, x0_id);
    x1 = netcdf.getVar(nc, x1_id);
    y0 = netcdf.getVar(nc, y0_id);
    y1 = netcdf.getVar(nc, y1_id);
    level = netcdf.getVar(nc,level_id);
    
  
    m1 = size(x1,1);
    m0 = size(x0,1);
    n1 = size(y1,1);
    n0 = size(y0,1);
    k  = size(level,1);
    
    fstagx = @(x) x(:);
    fstagy = @(y) y(1:(end-4));
    fgridx = @(x) x(2:(end-1));
    fgridy = @(y) y(2:(end-5));
    fgrid = @(M) M(2:(end-1),2:(end-5),:,:);
    %fstag = @(M) M(:,1:(end-3),:);
    fstag_layer = @(M) M(:,1:(end-4),:,:);
    
    data.xstag = fstagx(x0);
    data.xgrid = fgridx(x1);
    data.ystag = fstagy(y0);
    data.ygrid = fgridy(y1);
    data.level = level;
    data.time = time;
    
    data.thk = fgrid(netcdf.getVar(nc, thk_id,[0 0 0],[m1 n1 n_timeslices],[1 1 timestride]));
    data.thk_t = fgrid(netcdf.getVar(nc,thk_t_id,[0 0 0],[m1 n1 n_timeslices],[1 1 timestride]));
    data.bmlt = fgrid(netcdf.getVar(nc,bmlt_id,[0 0 0],[m1 n1 n_timeslices],[1 1 timestride]));
    data.lsurf = fgrid(netcdf.getVar(nc,lsurf_id,[0 0 0],[m1 n1 n_timeslices],[1 1 timestride]));
    if (have_flux_conv)
        uflx_conv = fgrid(netcdf.getVar(nc,uflx_conv_id,[0 0 0],[m1 n1 n_timeslices],[1 1 timestride]));
        vflx_conv = fgrid(netcdf.getVar(nc,vflx_conv_id,[0 0 0],[m1 n1 n_timeslices],[1 1 timestride]));
        data.flux_div = -uflx_conv-vflx_conv;
    end 
    
    data.uvel = fstag_layer(netcdf.getVar(nc, uvel_id,[0 0 0 0],[m0 n0 k n_timeslices],[1 1 1 timestride]));
    data.vvel = fstag_layer(netcdf.getVar(nc, vvel_id,[0 0 0 0],[m0 n0 k n_timeslices],[1 1 1 timestride]));
    
    data.uvelmean = squeeze(mean(data.uvel,3));
    data.vvelmean = squeeze(mean(data.vvel,3));
    
    %size(data.uvelmean)
    %size(data.stagthk)
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
      
  
    data.vel_div = zeros(size(data.thk,1),size(data.thk,2),size(data.thk,3));
    data.vel_div(2:(end-1),2:(end-1),:) = (ux(2:(end-1),2:(end-1),:) + ...
                                           vy(2:(end-1),2:(end-1),:)) .* ...
                                           data.thk(2:(end-1),2:(end-1),:);
    
    %data.flux_div = ...
    %     (data.uflx(2:end,    1:(end-1),:) + ...
    %      data.uflx(2:end,    2:end,:)- ...
    %      data.uflx(1:(end-1),1:(end-1),:)- ...
    %      data.uflx(1:(end-1),2:end,:))/(2*dx) + ...
    %     (data.vflx(1:(end-1),2:end,:)+ ...
    %      data.vflx(2:end,    2:end,:)- ...
    %      data.vflx(1:(end-1),1:(end-1),:)- ...
    %      data.vflx(2:end,    1:(end-1),:))/(2*dy);

    data.y_adv = zeros(size(data.thk,1),size(data.thk,2),size(data.thk,k));
    data.y_adv(:,2:(end-1),:) = ...
        data.vvelmean(round(m1/2),1,1) * ...
        (data.thk(:,3:end,:) - data.thk(:,1:(end-2),:))/(2*dy);
        
    [X1,Y1] = meshgrid(data.xgrid,data.ygrid);
    data.Xgrid = X1';
    data.Ygrid = Y1';
    

end
