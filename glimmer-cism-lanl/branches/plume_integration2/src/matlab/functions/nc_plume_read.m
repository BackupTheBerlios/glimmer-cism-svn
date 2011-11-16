function [data] = nc_plume_read(nc_filename,istart,timestride,iend)

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
    temp_id = netcdf.inqVarID(nc,'temp');
    salt_id = netcdf.inqVarID(nc,'salt');
    rhop_id = netcdf.inqVarID(nc,'rhop');
    
    x_id = netcdf.inqVarID(nc, 'x');
    y_id = netcdf.inqVarID(nc, 'y');

    time_dim_id = netcdf.inqDimID(nc,'time');
    time_var_id = netcdf.inqVarID(nc,'time');
    
    [tname,tlen] = netcdf.inqDim(nc,time_dim_id);
    if (iend < 0) 
      iend = tlen-1;
    end
 
    n_timeslices = floor((iend-istart)/timestride)+1;
    
    time = netcdf.getVar(nc,time_var_id,istart,n_timeslices,timestride);
    x = netcdf.getVar(nc, x_id);
    y = netcdf.getVar(nc, y_id);

    m = size(x,1);
    n = size(y, 1);

    bpos = netcdf.getVar(nc, bpos_id, [0 0 istart],[m n n_timeslices],[1 1 timestride]);
    pdep = netcdf.getVar(nc, pdep_id, [0 0 istart],[m n n_timeslices],[1 1 timestride]);
    bmelt = netcdf.getVar(nc, bmelt_id, [0 0 istart],[m n n_timeslices],[1 1 timestride]);
    train = netcdf.getVar(nc, train_id, [0 0 istart],[m n n_timeslices],[1 1 timestride]);
    entr = netcdf.getVar(nc, entr_id, [0 0 istart],[m n n_timeslices],[1 1 timestride]);
    su = netcdf.getVar(nc, su_id, [0 0 istart],[m n n_timeslices],[1 1 timestride]);
    sv = netcdf.getVar(nc, sv_id, [0 0 istart],[m n n_timeslices],[1 1 timestride]);
    u = netcdf.getVar(nc, u_id, [0 0 istart],[m n n_timeslices],[1 1 timestride]);
    v = netcdf.getVar(nc, v_id, [0 0 istart],[m n n_timeslices],[1 1 timestride]);
    temp = netcdf.getVar(nc, temp_id, [0 0 istart],[m n n_timeslices],[1 1 timestride]);
    salt = netcdf.getVar(nc, salt_id, [0 0 istart],[m n n_timeslices],[1 1 timestride]);
    rhop = netcdf.getVar(nc, rhop_id, [0 0 istart],[m n n_timeslices],[1 1 timestride]);
    
    netcdf.close(nc);

    waterdepth = max(max(bpos(:,:,1)));

    data.x = x((1+3+2):(end-3-2));
    data.y = y((1+3):(end-4));
    [xgrid,ygrid] = meshgrid(data.x,data.y);
    data.xgrid = xgrid';
    data.ygrid = ygrid';
    data.time = time;

    f = @(M) M((1+3+2):(end-3-2),(1+3):(end-4),:);

    data.bpos = f(bpos);
    data.pdep = f(pdep);
    data.bmelt = f(bmelt);
    data.train = f(train);
    data.entr = f(entr);
    data.su = f(su);
    data.sv = f(sv);
    data.u = f(u);
    data.v = f(v);
    data.temp = f(temp);
    data.salt = f(salt);
    data.rhop = f(rhop);
    
    data.draft = data.bpos - waterdepth;
    data.ambdepth = data.draft - data.pdep;
    [X,Y] = meshgrid(data.x,data.y);
    [data.gradx,data.grady,data.grad] = local_grad(X',Y',data.bpos);

end
