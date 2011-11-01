function [N] = pad_x(M)
  N = zeros(size(M,1)+2,size(M,2),size(M,3));
  N(2:end-1,:,:) =  M;
end

