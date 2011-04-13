function [data] = nc_plume_read(nc_filename,timeslice)

    nc = netcdf.open(nc_filename, 'NC_NOWRITE');
    
    su_id = netcdf.inqVarID(nc, 'su');
    sv_id = netcdf.inqVarID(nc, 'sv');
    u_id = netcdf.inqVarID(nc, 'u');
    v_id = netcdf.inqVarID(nc, 'v');
    time_id = netcdf.inqVarID(nc, 'time');
    bmelt_id = netcdf.inqVarID(nc, 'bmelt');
    bpos_id = netcdf.inqVarID(nc, 'bpos');    
    pdep_id = netcdf.inqVarID(nc, 'pdep');    
    train_id = netcdf.inqVarID(nc, 'train');
    entr_id = netcdf.inqVarID(nc, 'entr');

    x_id = netcdf.inqVarID(nc, 'x');
    y_id = netcdf.inqVarID(nc, 'y');
    
    time = netcdf.getVar(nc,time_id);
    x = netcdf.getVar(nc, x_id);
    y = netcdf.getVar(nc, y_id);

    bpos = netcdf.getVar(nc, bpos_id);
    pdep = netcdf.getVar(nc, pdep_id);
    bmelt = netcdf.getVar(nc, bmelt_id);
    train = netcdf.getVar(nc, train_id);
    entr = netcdf.getVar(nc, entr_id);
    su = netcdf.getVar(nc, su_id);
    sv = netcdf.getVar(nc, sv_id);
    u = netcdf.getVar(nc, u_id);
    v = netcdf.getVar(nc, v_id);
    
    if (timeslice < 0)
        timeslice = size(time,1);
    end
    waterdepth = max(max(bpos(:,:,1)));
    bmelt = bmelt(:,:,timeslice);
    bpos = bpos(:,:,timeslice);
    pdep = pdep(:,:,timeslice);
    train = train(:,:,timeslice);
    entr = entr(:,:,timeslice);
    su = su(:,:,timeslice);
    sv = sv(:,:,timeslice);
    u = u(:,:,timeslice);
    v = v(:,:,timeslice);

    data.x = x((1+3):(end-3));
    data.y = y((1+3):(end-4));
    f = @(M) M((1+3):(end-3),(1+3):(end-4),1);

    data.bpos = f(bpos);
    data.pdep = f(pdep);
    data.bmelt = f(bmelt);
    data.train = f(train);
    data.entr = f(entr);
    data.su = f(su);
    data.sv = f(sv);
    data.u = f(u);
    data.v = f(v);

    data.draft = data.bpos - waterdepth;
    data.ambdepth = data.draft - data.pdep;

end
