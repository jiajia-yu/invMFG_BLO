function [gradrhoL,gradmxL] = comp_gradL(siz,rho,mx,rho0,rho1,g,lambdaF,lambdaG)
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
dt = 1/siz(2);
logrho1 = log(rho1);
g = repmat(g,1,siz(2));

mxc = cat(1, mx(1,:)/2, (mx(1:end-1,:)+mx(2:end,:))/2, mx(end,:)/2);
msq = mxc.^2;

rhoc = cat(2, (rho0+rho(:,1))/2, (rho(:,1:end-1)+rho(:,2:end))/2);
logrho = log(rho);
rhoend  = rho(:,end);

% dynamic
ind = (rhoc>tol);

gradmxL = zeros(size(mxc));
gradmxL(ind) = dt*g(ind).*mxc(ind)./rhoc(ind);
gradmxL = (gradmxL(1:end-1,:)+gradmxL(2:end,:))/2;

gradrhoL = zeros(size(rhoc));
gradrhoL(ind) = -dt*g(ind).*msq(ind)./(2*rhoc(ind).^2);
gradrhoL = cat(2,(gradrhoL(:,1:end-1)+gradrhoL(:,2:end))/2,gradrhoL(:,end)/2);

% interaction
gradrhoL(:,1:end-1) = gradrhoL(:,1:end-1) + dt*lambdaF*(logrho(:,1:end-1)+1);

% terminal
ind = (rhoend>tol);
logterm = logrho(:,end)-logrho1;
gradrhoLend = gradrhoL(:,end);
gradrhoLend(ind) = gradrhoLend(ind) + lambdaG*(logterm(ind)+1);
gradrhoL(:,end) = gradrhoLend;

end