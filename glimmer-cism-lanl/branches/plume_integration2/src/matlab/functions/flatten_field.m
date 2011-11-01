function [flat_field] = flatten_field(gridded_field)

[m1,n1,k] = size(gridded_field);
m=m1-2;
n=n1-2;

remove_edge=@(d) d(2:(m1-1),2:(n1-1),:);
flat_field = reshape(remove_edge(gridded_field),m*n,1);

end
