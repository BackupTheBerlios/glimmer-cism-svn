function dfdy = Dy_xgrid(f,dy)

  dfdy = pad_edge(1,1,1/(4*dy)*(f(2:end,3:end  ,:)+f(1:end-1,3:end,:) ...
		     -1/(4*dy)*(f(2:end,1:end-2,:)+f(1:end-1,1:end-2,:))));

end
