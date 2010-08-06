
tauxy0 = [0,10,25,50];
acab = [-2.0, 0.0, 2.0];
js = ['1','2']

for k=1:length(js)
for i=1:length(tauxy0)
for j=1:length(acab)

filename = strcat('ssj',js(k),'_', sprintf('%i',tauxy0(i)),'_kPa_', ...
		  sprintf('%2.1f',acab(j)),'_acab');

[model_ythk, model_thk, model_yvel, model_vvel] = read_steadyice( sprintf('%s%s',filename,'.out.nc')  );

save( sprintf('%s%s',filename,'_model.mat'), 'model_ythk','model_thk','model_yvel','model_vvel');

end 
end
end



