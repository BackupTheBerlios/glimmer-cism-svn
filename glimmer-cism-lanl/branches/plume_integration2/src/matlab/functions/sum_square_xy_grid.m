function ans=sum_square_xy_grid(u,v)

  ans = pad_edge(1,1,...
      0.25*(u(1:end-1,2:end-1,:)+u(2:end,2:end-1,:)).^2 + ...
      0.25*(v(2:end-1,1:end-1,:)+v(2:end-1,2:end,:)).^2 );

end
