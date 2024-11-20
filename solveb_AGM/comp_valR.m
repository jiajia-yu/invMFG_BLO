function valR = comp_valR(b)
nx = size(b,1);
ny = size(b,2);

grad1 = (b(2:end,:)-b(1:end-1,:));
grad1 = (grad1(:,1:end-1)+grad1(:,2:end))/2;
grad2 = (b(:,2:end)-b(:,1:end-1));
grad2 = (grad2(1:end-1,:)+grad2(2:end,:))/2;
valR = 0.5*sum(grad1.^2*(nx/ny) + grad2.^2*(ny/nx),'all');

end