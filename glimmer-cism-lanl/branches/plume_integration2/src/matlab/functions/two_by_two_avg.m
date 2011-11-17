function res = two_by_two_avg(M)

  [m,n] = size(M);

res = zeros(m,n);

for i=1:(m-1)
  for j=1:(n-1)
    res(i,j) = (1/4)*(sum(sum(M(i:i+1,j:j+1))));
end
end


end
