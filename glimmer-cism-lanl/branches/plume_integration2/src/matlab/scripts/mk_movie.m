

fid = fopen('/home/gehne/Work/Programming/Spectra_Z/era-40_2.5-deg_Z-129_200-hPa_6hr-anl_19830701-19830710.dat');
%fid = fopen('/home/gehne/Work/Programming/Spectra_Z/era-40_2.5-deg_Z-129_925-hPa_6hr-anl_198307-200208.dat');
D10days = fread(fid,[144,2920],'single');
fclose(fid);

M1 = zeros(40,17,144);

for i=1:40
        M1(i,:,:) = D10days(:,45+(i-1)*73:-1:29+(i-1)*73)';
end

%% find the reconstruction at each time
nPCF=20;
P = ParabolCyl([-100:2.5:100]',nPCF-1);
P20 = ParabolCyl([-20:2.5:20]',nPCF-1);

[m,n]=size(P);




%% make the movie

lat=[-20:2.5:20];
lon=[0:2.5:357.5];
h = figure('Position',[300 300 800 300]);

subplot(3,1,1);
contourf(lon,lat,squeeze(M1(1,:,:)),10,'EdgeColor','none');
%contourf(lon,lat,squeeze(M1(1,:,:)));
colormap jet;
colorbar;
subplot(3,1,2);
subplot(3,1,3);
caxis(caxis) 


reruns=1;                  % number of times movie is to play
fps=4;                     % frames per second

nframes =40;              % number of frames in the movie
Frames = moviein(nframes); % initialize the matrix 'Frames'

    

for j = 1:nframes
    data = squeeze(M1(j,:,:));
    mean_dat = sum(sum(data,1),2)/(17*144);


%% pad data with b value outside of +-20 degrees
%    b = mean_dat;
    b = 121000;
    data_pad = b*ones(m,144);

%data_pad = zeros(81,144);
    data_pad(floor(m/2)-7:floor(m/2)+9,:) = data;

%calculate the mean in each longitude and subtract
%    mean = sum(data_pad,1)/m;
%    one = ones(m,1);

%    data_pad= data_pad-one*mean;

% subtract total mean
    data_pad=data_pad-b;

%% solve the linear system P'P*a = P'*data_pad to find the coefficients a

    a = (P'*P)\(P'*data_pad);

%% reconstruct the data

    rec = P20*a;
%    rec=rec+ones(17,1)*mean; 
    rec=rec+b;
    
%% find the RMSE between data and reconstruction

    RMSE = sqrt(sum(abs(data-rec).^2,1)/17);

    maxerr = max(RMSE);
    
    disp(['the mean of the data at this time is ' num2str(mean_dat)]);
    
    subplot(3,1,1,'replace');
    contourf(lon,lat,squeeze(M1(j,:,:)),10,'EdgeColor','none');
    title(['Geopotential 850mb -- day =',num2str(j/4)])
    colormap jet;
    colorbar;
    subplot(3,1,2,'replace');
    contourf(lon,lat,rec,10,'EdgeColor','none');
    title(['Geopotential 850mb reconstruction, b = ', num2str(b)]);
    colormap jet;
    colorbar;
    subplot(3,1,3,'replace')
    plot([0:2.5:357.5],RMSE);
    xlabel('longitude');
    ylabel('error');
    title('RMSE');

    
    
   
    Frames(j) = getframe(h);

end   

 
 map=colormap;    % Uses the previously defined colormap 
 movie2avi(Frames,'test','fps',fps)