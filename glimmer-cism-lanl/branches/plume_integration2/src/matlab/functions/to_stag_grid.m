function N=to_stag_grid(M)
  N = M(2:end,2:end,:)+M(1:end-1,2:end,:)+M(2:end,1:end-1,:)+M(1:end-1,1:end-1,:);
end

