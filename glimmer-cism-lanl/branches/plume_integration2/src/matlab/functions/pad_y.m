function [N] = pad_y(M)
  N = zeros(size(M,1),size(M,2)+2,size(M,3));
  N(:,2:end-1,:) = M;

end
