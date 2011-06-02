function [ ug,vg, uag, vag, agfrac ] = geostrophic_vel( d_plume )

    g = 9.81;
    betaS = 7.86e-4;
    betaT = 3.87e-5;
    f = 2*(2*pi)/(3600.0*24.0)*sin(80 * pi/180);
    
    Sa = 34.75;
    Ta = 0.0;
    
    dx = d_plume.x(2)-d_plume.x(1);
    dy = d_plume.y(2)-d_plume.y(1);
    
    salt = d_plume.salt;
    temp = d_plume.temp;
    ambsurf = d_plume.bpos - d_plume.pdep;
    pdep = d_plume.pdep;
    
    [m,n,k] = size(d_plume.su);
    
    ug = zeros(m,n,k);
    vg = zeros(m,n,k);
    
    agfrac = zeros(m,n,k);
    
    ug(:,2:(end-1),:) = (g/f)*( ...
          0.5*pdep(:,2:(end-1),:) .* ...
               (betaS*(salt(:,3:end,:)-salt(:,1:(end-2),:))/(2*dy) ...
              - betaT*(temp(:,3:end,:)-temp(:,1:(end-2),:))/(2*dy)) ...
            +  (betaS*(Sa - salt(:,2:(end-1),:))-betaT*(Ta-temp(:,2:(end-1),:))) ...
             .*(ambsurf(:,3:end,:)-ambsurf(:,1:(end-2),:))/(2*dy));
    
    vg(2:(end-1),:,:) = (-g/f)*( ...
          0.5*pdep(2:(end-1),:,:) .* ...
               (betaS*(salt(3:end,:,:)-salt(1:(end-2),:,:))/(2*dx) ...
              - betaT*(temp(3:end,:,:)-temp(1:(end-2),:,:))/(2*dx)) ...
             +  (betaS*(Sa - salt(2:(end-1),:,:))-betaT*(Ta-temp(2:(end-1),:,:))) ...
             .*(ambsurf(3:end,:,:)-ambsurf(1:(end-2),:,:))/(2*dx));
         
     
    
    uag = zeros(m,n,k);
    vag = zeros(m,n,k);
    
    uag(:,2:(end-1),:) = d_plume.su(:,2:(end-1),:) - ug(:,2:(end-1),:);
    vag(2:(end-1),:,:) = d_plume.sv(2:(end-1),:,:) - vg(2:(end-1),:,:);
   
    cutedge = @(M) M(2:(end-1),2:(end-1),:);
    
    agfrac(2:(end-1),2:(end-1),:) = ...
          (cutedge(uag).^2 + cutedge(vag).^2) ./  ...
          (cutedge(ug).^2 + cutedge(vg).^2);
           
    
                
end

