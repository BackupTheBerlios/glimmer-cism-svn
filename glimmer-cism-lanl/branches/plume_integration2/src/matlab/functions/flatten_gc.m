function [ flat_ocean, flat_ice ] = flatten_gc( docean,dice, times)

[m1,n1,k] = size(docean.bmelt(:,:,times));

k2 = length(dice.time(times));
if (k ~= k2)
    error('number of timeslices disagree in ice and plume data')
end

m=m1-2;
n=n1-2;

remove_edge=@(d) d(2:(m1-1),2:(n1-1),:);
%remove_edge=@(d) d;

%extra_rows = 40;
%remove_hot_spot=@(d) d(2:(m1-1),(2+extra_rows):(n1-1),:);

ygridrep = repmat(docean.ygrid,[1,1,k]);

secondderivDxDx = zeros(m1,n1,k);

dx = docean.xgrid(2) - docean.xgrid(1);

for ti=1:k

 for xi=3:(m1-2)
  for yi=1:n1
  p = polyfit(dx*[-2,-1,0,1,2],docean.draft((xi-2):(xi+2),yi,ti)',2);
  secondderivDxDx(xi,yi,ti) = p(1);
  end
end
end

%secondderivDxDx(2:(m1-1),:,times) = docean.draft(3:m1,:,times)+ ...
%                                   docean.draft(1:(m1-2),:,times) - ...
%
                                 2*docean.draft(2:(m1-1),:,times);

flat_ocean.curvature = reshape(remove_edge(secondderivDxDx(:,:,1:k)),m*n*k,1);
flat_ocean.y    =      reshape(remove_edge(ygridrep),m*n*k,1);
flat_ocean.bmlt = reshape(remove_edge(docean.bmelt(:,:,times)),m*n*k,1);
flat_ocean.u =    reshape(remove_edge(docean.u(:,:,times)),m*n*k,1);
flat_ocean.v =    reshape(remove_edge(docean.v(:,:,times)),m*n*k,1);
flat_ocean.su =   reshape(remove_edge(docean.su(:,:,times)),m*n*k,1);
flat_ocean.sv =   reshape(remove_edge(docean.sv(:,:,times)),m*n*k,1);
flat_ocean.speed = reshape(sqrt(remove_edge(docean.su(:,:,times)) .^2 ...
                               +remove_edge(docean.sv(:,:,times)) .^ 2),m*n*k,1);
flat_ocean.temp = reshape(remove_edge(docean.temp(:,:,times)),m*n*k,1);
flat_ocean.salt = reshape(remove_edge(docean.salt(:,:,times)),m*n*k,1);
flat_ocean.train =reshape(remove_edge(docean.train(:,:,times)/1000.0),m*n*k,1);
flat_ocean.gradx = reshape(remove_edge(docean.gradx(:,:,times)),m*n*k,1);
flat_ocean.grady = reshape(remove_edge(docean.grady(:,:,times)),m*n*k,1);
flat_ocean.grad = reshape(remove_edge(docean.grad(:,:,times)),m*n*k,1);
flat_ocean.draft = reshape(remove_edge(docean.draft(:,:,times)),m*n*k,1);

flat_ice.y_adv = reshape(remove_edge(-dice.y_adv(:,:,times)),m*n*k,1);                           
flat_ice.ice_def = reshape(remove_edge(-dice.flux_div(:,:,times)+dice.y_adv(:,:,times)),m*n*k,1);
flat_ice.thk_div = reshape(remove_edge(dice.vel_div(:,:,times)),m*n*k,1);
flat_ice.flux_div = reshape(remove_edge(dice.flux_div(:,:,times)),m*n*k,1);

%flat_ice.downstream_ice_def = reshape(remove_hot_spot(-dice.flux_div(:,:,times)+dice.y_adv(:,:,times)), ...
%                                      m*(n-extra_rows)*k,1);
%flat_ice.downstream_draft = reshape(remove_hot_spot(docean.draft(:,:,times)), ...
%                                    m*(n-extra_rows)*k,1);


[m3,n3,k3] = size(docean.draft);
channel_amp = zeros(m3,n3,k3);

for i=1:n3
  channel_amp(:,i,:) = docean.draft(:,i,:) - repmat(min(docean.draft(:,i,:)),[m3 1]);
end

%flat_ice.downstream_channel_amp = reshape(remove_hot_spot(channel_amp(:,:,times)),...
%                                    m*(n-extra_rows)*k,1);

end

