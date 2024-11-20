function [valL,gradrhoL,gradmxL,gradmyL] = comp_valgradL(siz,rho,mx,my,rho0,rho1,b,lambdaF,lambdaG)
% [valL,gradrhoL,gradmL] = comp_valgradL(siz,rho,m,rho0)
% compute the value of lower level objective L
% and the gradient of lower level objective L
% with respect to density rho and flux m
% siz = [nx,ny,nt]
% rho, density, matrix of size (nx,ny,nt)
% m, flux, cell of size 1*2,
%   m{1}, matrix of size (nx-1,ny,nt)
%   m{2}, matrix of size (nx,ny-1,nt)
% rho0, initial density, matrix of size (nx,ny)
% rho1, desired terminal density, matrix of size (nx,ny)
% b, interaction cost, matrix of size (nx,ny)
% the values are normalized by temporal grid size

% precomputation
tol = 1e-8;
dt = 1/siz(3);
logrho1 = log(rho1);

mxc = cat(1, mx(1,:,:)/2, (mx(1:end-1,:,:)+mx(2:end,:,:))/2, mx(end,:,:)/2);
myc = cat(2, my(:,1,:)/2, (my(:,1:end-1,:)+my(:,2:end,:))/2, my(:,end,:)/2);
msq = mxc.^2 + myc.^2;

rhoc = cat(3, (rho0+rho(:,:,1))/2, (rho(:,:,1:end-1)+rho(:,:,2:end))/2);
logrho = log(rho);
rhoend  = rho(:,:,end);

% dynamic
ind = (rhoc>tol);
valL = 0.5*sum(msq(ind)./rhoc(ind));

gradmxL = zeros(size(mxc));
gradmxL(ind) = dt*mxc(ind)./rhoc(ind);
gradmxL = (gradmxL(1:end-1,:,:)+gradmxL(2:end,:,:))/2;
gradmyL = zeros(size(myc));
gradmyL(ind) = dt*myc(ind)./rhoc(ind);
gradmyL = (gradmyL(:,1:end-1,:)+gradmyL(:,2:end,:))/2;

gradrhoL = zeros(size(rhoc));
gradrhoL(ind) = -dt*msq(ind)./(2*rhoc(ind).^2);
gradrhoL = cat(3,(gradrhoL(:,:,1:end-1)+gradrhoL(:,:,2:end))/2,gradrhoL(:,:,end)/2);

% interaction
valL = dt*(valL + sum(rho(:,:,1:end-1).*b,'all') ...
        + lambdaF*sum(rho(:,:,1:end-1).*logrho(:,:,1:end-1),'all'));

gradrhoL(:,:,1:end-1) = gradrhoL(:,:,1:end-1) + dt*b + dt*lambdaF*(logrho(:,:,1:end-1)+1);

% terminal
ind = (rhoend>tol);
logterm = logrho(:,:,end)-logrho1;
valL = valL + lambdaG*sum(rhoend(ind).*logterm(ind));
gradrhoLend = gradrhoL(:,:,end);
gradrhoLend(ind) = gradrhoLend(ind) + lambdaG*(logterm(ind)+1);
gradrhoL(:,:,end) = gradrhoLend;


end