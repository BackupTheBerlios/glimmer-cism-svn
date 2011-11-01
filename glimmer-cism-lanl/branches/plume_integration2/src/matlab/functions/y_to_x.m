function N = y_to_x(M)

[m,n,k] = size(M);
N = zeros(m-1,n+1,k);
N(1:(m-1),1  ,:) = 0.5*(M(1:(m-1),1,:)+M((2:m),1,:));
N(1:(m-1),n+1,:) = 0.5*(M(1:(m-1),n,:)+M((2:m),n,:));
N(1:(m-1),2:n,:) = 0.25*(M(1:(m-1),2:n,:)+M(2:m,2:n,:)+M(1:(m-1),1:(n-1),:)+M(2:m,1:(n-1),:));
  
end
