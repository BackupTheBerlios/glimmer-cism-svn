%data is in variable d

%f = '/data/gladish/may3_coupled_super_channels_trial_90/plume.may3_coupled_super_channels_trial_90.out.nc';
f_may3 = '/data/gladish/may3.1_diff_10.0_cdb_2.5_dtime_14400.0_entype_6_pthk_10.0_tvel_0.0/plume.may3.1_diff_10.0_cdb_2.5_dtime_14400.0_entype_6_pthk_10.0_tvel_0.0.out.nc';

d3 = nc_plume_read(f_may3);

fs = 12;

fig1 = figure('Position',[300 300 600 400]);
pause(1);
winsize = get(fig1,'Position');
m = size(d3.su,1);
n = size(d3.su,2);
numframes=size(d.su,3);
numframes = 330;
%mframes=moviein(numframes,fig1,winsize);

set(fig1,'NextPlot','replacechildren');

channel_amp = zeros(m,n);
mean_draft = zeros(m,n);
max_draft = zeros(m,n);

x = d3.x;
y = d3.y;

for k=1:numframes
    
    sprintf('%d %d',k,numframes)
    
    draft = squeeze(d3.draft(:,:,k));
    bpos = squeeze(d3.bpos(:,:,k));
    
    for i=1:n
        channel_amp(:,i) = draft(:,i) - min(draft(:,i));
        mean_draft(:,i) = sum(draft(:,i))/m;
        max_draft(:,i) = min(draft(:,i));
    end
    
    %hold on
    subplot(1,2,2);
    contourf(x/1000.0,y/1000.0,channel_amp',30,'EdgeColor','None') ;colorbar;
    %[C,h] = contour(x/1000.0,y/1000.0,max_draft',20,'w');
    %caxis([min(min(channel_amp)) max(max(channel_amp))]);
    caxis([0,300.0]);
    xlabel('Across shelf distance (km)','FontSize',fs);
    ylabel('Along shelf distance (km)','FontSize',fs);
    %clabel(C,'FontSize',fs,'Color','k')
    %title('Channel depth with deepest-draft contours','FontSize',fs);
    title('Channel depth','FontSize',fs);
    
    subplot(1,2,1);
    contourf(x/1000.0,y/1000.0,draft',30,'EdgeColor','None');colorbar;
    caxis([-600,0.0]);
    xlabel('Across shelf distance (km)','FontSize',fs);
    ylabel('Along shelf distance (km)','FontSize',fs);
    %clabel(C,'FontSize',fs,'Color','k')
    %title('Channel depth with deepest-draft contours','FontSize',fs);
    title('Ice draft','FontSize',fs);
    
    %A(:,i)=getframe(fig1,winsize);
    mframes(k) = getframe(gcf);
    
end
