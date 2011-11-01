function [data_avg] = nc_plume_avg(data)

data_avg.x = data.x;
data_avg.y = data.y;
data_avg.xgrid = data.xgrid;
data_avg.ygrid = data.ygrid;
data_avg.time = [data.time(end)];
data_avg.bpos = mean(data.bpos,3);
data_avg.pdep = mean(data.pdep,3);
data_avg.bmelt = mean(data.bmelt,3);
data_avg.train = mean(data.train,3);
data_avg.entr = mean(data.entr,3);
data_avg.su = mean(data.su,3);
data_avg.sv = mean(data.sv,3);
data_avg.u = mean(data.u,3);
data_avg.v = mean(data.v,3);
data_avg.temp = mean(data.temp,3);
data_avg.salt = mean(data.salt,3);
data_avg.rhop = mean(data.rhop,3);
data_avg.draft = mean(data.draft,3);
data_avg.ambdepth = mean(data.draft-data.pdep,3);
data_avg.gradx = mean(data.gradx,3);
data_avg.grady = mean(data.grady,3);
data_avg.grad = mean(data.grad,3);

end 
