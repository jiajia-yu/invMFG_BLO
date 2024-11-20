function gradbR = comp_gradR(b)

nx = size(b,1);
ny = size(b,2);

grad1 = (b(2:end,:)-b(1:end-1,:));
grad1 = (grad1(:,1:end-1)+grad1(:,2:end))/2;
grad2 = (b(:,2:end)-b(:,1:end-1));
grad2 = (grad2(1:end-1,:)+grad2(2:end,:))/2;

gradbR1 = conv2(grad1,[-1,-1;1,1],'full') * (0.5*nx/ny);
gradbR2 = conv2(grad2,[-1,1;-1,1],'full') * (0.5*ny/nx);
% gradbR1 = zeros(nx,ny);
% gradbR1(1:end-1,1:end-1) = - grad1;
% gradbR1(1:end-1,2:end) = gradbR1(1:end-1,2:end) - grad1;
% gradbR1(2:end,1:end-1) = gradbR1(2:end,1:end-1) + grad1;
% gradbR1(2:end,2:end) = gradbR1(2:end,2:end) + grad1;
% gradbR1 = gradbR1 * (0.5*nx/ny);

gradbR = gradbR1 + gradbR2;


end