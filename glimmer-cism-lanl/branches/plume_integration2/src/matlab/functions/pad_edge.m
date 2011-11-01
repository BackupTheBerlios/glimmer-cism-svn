function [N] = pad_edge(m,n,M)
% pad the first dimension with m cells on either side
% pad the second dimension with n cells on either side
N = [zeros(m,size(M,2)+2*n,size(M,3)); ...
     zeros(size(M,1),n,size(M,3)), M,  zeros(size(M,1),n,size(M,3)); ...
     zeros(m,size(M,2)+2*n,size(M,3))];
end
