function dfdx = Dx_ygrid(f)

  dfdx = pad_edge(1,1,0.5*(f(3:end  ,2:end,:)+f(3:end  ,1:end-1,:) ...
	  	     -0.5*(f(1:end-2,2:end,:)+f(1:end-2,1:end-1,:))));

end
