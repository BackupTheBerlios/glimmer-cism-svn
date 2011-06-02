function [gradx,grady,grad] = local_grad(x,y,phi)

[m,n,k] = size(phi);

gradx = zeros(m,n,k);
grady = zeros(m,n,k);
grad = zeros(m,n,k);

for i=2:(m-1)
    for j=2:(n-1)
        gradx(i,j,:) = (phi(i+1,j,:)-phi(i-1,j,:))/(x(i+1,j)-x(i-1,j));
        grady(i,j,:) = (phi(i,j+1,:)-phi(i,j-1,:))/(y(i,j+1)-y(i,j-1));
        grad(i,j,:) = sqrt(gradx(i,j,:).^2 + grady(i,j,:).^2);
    end
end
end