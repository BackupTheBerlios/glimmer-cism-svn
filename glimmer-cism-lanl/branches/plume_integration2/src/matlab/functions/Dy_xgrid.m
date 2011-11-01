function dfdy = Dy_xgrid(f)

  dfdy = pad_edge(1,1,0.5*(f(2:end,3:end  ,:)+f(1:end-1,3:end,:) ...
		     -0.5*(f(2:end,1:end-2,:)+f(1:end-1,1:end-2,:))));

end
