function dfdx = Dx_ygrid(f,dx)

  dfdx = pad_edge(1,1,  1/(4*dx)*(f(3:end  ,2:end,:)+f(3:end  ,1:end-1,:) ...
    	               -1/(4*dx)*(f(1:end-2,2:end,:)+f(1:end-2,1:end-1,:))));

end
