function [dfdx] = Dx_xgrid(f)
  dfdx = pad_edge(1,0, f(2:end,:,:)-f(1:end-1,:,:));
end
