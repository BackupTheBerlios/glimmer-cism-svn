function res = three_by_three_avg(M)

  [m,n] = size(M);

res = zeros(m,n);

for i=2:(m-1)
  for j=2:(n-1)
    res(i,j) = (1/9)*(sum(sum(M(i-1:i+1,j-1:j+1))));
end
end


end
