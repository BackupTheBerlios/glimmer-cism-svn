function [data] = nc_ice_read(nc_filename, istart,timestride, iend,minslices)

    nc = netcdf.open(nc_filename, 'NC_NOWRITE');
       
    thk_id = netcdf.inqVarID(nc, 'thk');
    thk_t_id = netcdf.inqVarID(nc,'thk_t');

    efvs_id =  netcdf.inqVarID(nc,'efvs');
    lsurf_id = netcdf.inqVarID(nc,'lsurf');
    uvel_id = netcdf.inqVarID(nc, 'uvelhom');
    vvel_id = netcdf.inqVarID(nc, 'vvelhom');
    bmlt_id = netcdf.inqVarID(nc, 'bmlt');
    x0_id = netcdf.inqVarID(nc, 'x0');
    y0_id = netcdf.inqVarID(nc, 'y0');
    x1_id = netcdf.inqVarID(nc, 'x1');
    y1_id = netcdf.inqVarID(nc, 'y1');
    level_id =netcdf.inqVarID(nc,'level');
    time_dim_id = netcdf.inqDimID(nc,'time');
    time_var_id = netcdf.inqVarID(nc,'time');

    try
    H_diff_id = netcdf.inqVarID(nc,'H_diff_t');
        have_H_diff = true;
    catch
        have_H_diff = false;
    end

    try
        uflx_conv_id = netcdf.inqVarID(nc,'uflx_conv');
        vflx_conv_id = netcdf.inqVarID(nc,'vflx_conv');
        have_flux_conv = true;
    catch 
        have_flux_conv = false;
    end
    
    [tname,tlen] = netcdf.inqDim(nc,time_dim_id);
    if (istart < 0)
      if (iend > 0) 
	     error('Does not make send to have istart < 0 and iend > 0');
      end
      istart = tlen-minslices;
      stride = 1;
    end
    if (iend < 0)
      iend = tlen-1;
    end
   
    n_timeslices = floor((iend-istart)/timestride)+1;

    time = netcdf.getVar(nc,time_var_id,istart,n_timeslices,timestride);
    
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
    
    fstagx = @(x) x(3:m0-2);
    fstagy = @(y) y(1:n0-3);
    fstag_layer = @(M) M(3:m0-2,1:n0-3,:,:);

    fgridx = @(x) x(4:m1-3);
    fgridy = @(y) y(2:n1-4);
    fgrid = @(M) M(4:m1-3,2:n1-4,:,:);

    
    data.xstag = fstagx(x0);
    data.xgrid = fgridx(x1);
    data.ystag = fstagy(y0);
    data.ygrid = fgridy(y1);
    data.level = level;
    data.time = time;
    
    data.thk = fgrid(netcdf.getVar(nc, thk_id,[0 0 istart],[m1 n1 n_timeslices],[1 1 timestride]));
    data.thk_t = fgrid(netcdf.getVar(nc,thk_t_id,[0 0 istart],[m1 n1 n_timeslices],[1 1 timestride]));

if (have_H_diff)
    data.H_diff_t = fgrid(netcdf.getVar(nc,H_diff_id,[0 0 istart],[m1 n1 n_timeslices],[1 1 timestride]));
end

    data.bmlt = fgrid(netcdf.getVar(nc,bmlt_id,[0 0 istart],[m1 n1 n_timeslices],[1 1 timestride]));
    data.lsurf = fgrid(netcdf.getVar(nc,lsurf_id,[0 0 istart],[m1 n1 n_timeslices],[1 1 timestride]));
    data.efvs = fgrid(netcdf.getVar(nc,efvs_id,[0 0 0 istart],[m1 n1 k n_timeslices],[1 1 1 timestride]));

    if (have_flux_conv)
        uflx_conv = fgrid(netcdf.getVar(nc,uflx_conv_id,[0 0 istart],[m1 n1 n_timeslices],[1 1 timestride]));
        vflx_conv = fgrid(netcdf.getVar(nc,vflx_conv_id,[0 0 istart],[m1 n1 n_timeslices],[1 1 timestride]));
        data.flux_div = -uflx_conv-vflx_conv;
data.uflx_conv = uflx_conv;
data.vflx_conv = vflx_conv;
    end 
    
    data.uvel = fstag_layer(netcdf.getVar(nc, uvel_id,[0 0 0 istart],[m0 n0 k n_timeslices],[1 1 1 timestride]));
    data.vvel = fstag_layer(netcdf.getVar(nc, vvel_id,[0 0 0 istart],[m0 n0 k n_timeslices],[1 1 1 timestride]));
    
netcdf.close(nc);

    data.uvelmean = squeeze(mean(data.uvel,3));
    data.vvelmean = squeeze(mean(data.vvel,3));
    
    
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

      size(data.uvelmean);
uy = (data.uvelmean(1:(end-1),2:end    ,:) + ...
      data.uvelmean(2:end    ,2:end    ,:) - ...
      data.uvelmean(1:(end-1),1:(end-1),:) - ...
      data.uvelmean(2:end    ,1:(end-1),:))/(2*dx);

vx = (data.vvelmean(2:end    ,1:(end-1),:)+ ...
      data.vvelmean(2:end    ,2:end    ,:)- ...
      data.vvelmean(1:(end-1),1:(end-1),:)- ...
      data.vvelmean(1:(end-1),2:end    ,:))/(2*dx);

    data.vel_div = zeros(size(data.thk,1),size(data.thk,2),size(data.thk,3));
    data.vel_div(2:(end-1),2:(end-1),:) = (ux(2:(end-1),2:(end-1),:) + ...
                                           vy(2:(end-1),2:(end-1),:)) .* ...
                                           data.thk(2:(end-1),2:(end-1),:);
    
data.vorticity = zeros(size(data.thk,1),size(data.thk,2),size(data.thk,3));
size(uy);
size(vx);
size(data.vorticity);
size(data.thk);
data.vorticity = vx - uy;

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
