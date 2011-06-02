function [x0,y0,x1,y1,thk,uvel,vvel] = nc_read(nc_filename,timeslice)

    nc = netcdf.open(nc_filename, 'NC_NOWRITE');
    
    thk_id = netcdf.inqVarID(nc, 'thk');
    uvel_id = netcdf.inqVarID(nc, 'uvelhom');
    vvel_id = netcdf.inqVarID(nc, 'vvelhom');
    
    x0_id = netcdf.inqVarID(nc, 'x0');
    y0_id = netcdf.inqVarID(nc, 'y0');
    x1_id = netcdf.inqVarID(nc, 'x1');
    y1_id = netcdf.inqVarID(nc, 'y1');
    
    x0 = netcdf.getVar(nc, x0_id);
    x1 = netcdf.getVar(nc, x1_id);
    y0 = netcdf.getVar(nc, y0_id);
    y1 = netcdf.getVar(nc, y1_id);
    thk = netcdf.getVar(nc, thk_id);
    uvel = netcdf.getVar(nc, uvel_id);
    vvel= netcdf.getVar(nc, vvel_id);
    
    thk = thk(:,:,timeslice);
    uvel = uvel(:,:,:,timeslice);
    vvel = vvel(:,:,:,timeslice);
    
end