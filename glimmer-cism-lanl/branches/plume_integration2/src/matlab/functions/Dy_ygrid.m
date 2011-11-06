function [dfdy] = Dy_ygrid(f,dy)
  dfdy = pad_edge(0,1, (f(:,2:end,:) -f(:,1:end-1,:))/dy);
end
