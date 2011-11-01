function [data_avg] = nc_ice_avg(data)

data_avg.xstag = data.xstag;
data_avg.ystag = data.ystag;
data_avg.xgrid = data.xgrid;
data_avg.ygrid = data.ygrid;
data_avg.level = data.level;
data_avg.time  = [data.time(end)];

data_avg.thk = mean(data.thk,3);
data_avg.thk_t = mean(data.thk_t,3);
data_avg.bmlt = mean(data.bmlt,3);
data_avg.lsurf =mean(data.lsurf,3);
data_avg.efvs = mean(data.efvs,4);
%data_avg.uflx_conv = mean(data.uflx_conv,4);
%data_avg.vflx_conv = mean(data.vflx_conv,4);
data_avg.flux_div = mean(data.flux_div,4);
data_avg.uvel = mean(data.uvel,4);
data_avg.vvel = mean(data.vvel,4);
data_avg.uvelmean = mean(data.uvelmean,3);
data_avg.vvelmean = mean(data.vvelmean,3);
data_avg.vel_div = mean(data.vel_div,3);
data_avg.y_adv = mean(data.y_adv,3);
data_avg.Xgrid = data.Xgrid;
data_avg.Ygrid = data.Ygrid;

end
