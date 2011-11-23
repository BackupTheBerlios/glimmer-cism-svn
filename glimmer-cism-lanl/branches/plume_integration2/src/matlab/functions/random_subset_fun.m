function fun  = random_subset_fun(f,N)

%%%% we are going to return a lambda-style function that pulls out a
%%%% random subset of a given vector v

N2 = floor(f*N);
indices = zeros(N2,1);

for i=1:N2
  found_new = false;
  while (not(found_new))
    irand = floor(rand*N);
    found_new = true;
    for j=1:(i-1)
      if (indices(j) == irand)
	found_new = false;
      end
    end
  end
    indices(i) = irand;
end

indices = min(indices+1,N);

fun = @(v) v(indices);


end

