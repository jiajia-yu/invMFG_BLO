function [rho,mx,my,projerr] = proj_rhom(rho,mx,my,rho0)
% [rho,m] = proj_rhom(rho,m)
% [rho,m] = proj_rhom(rho,m,rho0)
% solves the projection appears in 2-D potential mean-field games
% rho: density, matrix of size nx*ny*nt
% m: flux, cell of size 1*2
%         m{1}, matrix of size (nx-1)*ny*nt
%         m{2}, matrix of size nx*(ny-1)*nt
% rho0: initial density, matrix of size nx*ny
%       if not given, treat as zeros(nx,ny)

persistent siz sizsc lapinv

% check if the variable size has changed
% and we need to precompute eigenvalues
precompute=0;
if  ~exist('siz','var')||numel(siz)==0
    precompute=1;
elseif sum(abs(size(rho)-siz))>0
    precompute=1;
end

% precompute eigenvalues for inverting negative laplacian operator
if precompute
    % size
    siz = size(rho);
    % scale to reduce operator norm
    sizmax = max(siz);
    sizsc = siz/sizmax;
    % eigenvalues
    lap_x = reshape((2 - 2*cos(pi*(0:siz(1)-1)/siz(1)))...
        *(sizsc(1)*sizsc(1)),siz(1),1,1);
    lap_y = reshape((2 - 2*cos(pi*(0:siz(2)-1)/siz(2)))...
        *(sizsc(2)*sizsc(2)),1,siz(2),1);
    lap_t = reshape((2 - 2*cos(pi*(1:2:2*siz(3)-1)/(2*siz(3)+1)))...
        *(sizsc(3)*sizsc(3)),1,1,siz(3));
    lap = repmat(lap_x,1,siz(2),siz(3)) ...
        + repmat(lap_y,siz(1),1,siz(3)) ...
        + repmat(lap_t,siz(1),siz(2),1);
    lapinv = 1./lap;
end

if (nargin < 4)
    rho0 = zeros(siz(1:2));
elseif sum(abs(size(rho0)-siz(1:2)))>0
    error('Inconsistent input');
end

% %-----check op_dt and op_dtadj
% dtrho = op_dt(rho);
% phi = rand(size(dtrho));
% dtadjphi = op_dtadj(phi);
% err = sum(rho.*dtadjphi,'all')-sum(dtrho.*phi,'all')
% %-----check op_div and op_divadj
% divm = op_div(m);
% phi = rand(size(divm));
% divadjphi = op_divadj(phi);
% divmphi = sum(divm.*phi,'all');
% mdivadjphi = sum(m{1}.*divadjphi{1},'all')+sum(m{2}.*divadjphi{2},'all');
% err = divmphi - mdivadjphi

phi = op_dt(rho) + op_div(mx,my);
phi(:,:,1) = phi(:,:,1) - sizsc(1)*rho0;
phi = lapinv.*lsbd_dct(phi);
phi = lsbd_idct(phi);

rho = rho - op_dtadj(phi);
[divadjphix,divadjphiy] = op_divadj(phi);
mx = mx - divadjphix;
my = my - divadjphiy;
% m = cellfun(@minus,m,op_divadj(phi),'UniformOutput',false);

%-----check whether op_dt(rho) + op_div(m) 
%                  -sizsc(1)*cat(1,rho0,zeros(nt-1,nx+1,ny+1))= 0
projerr = op_dt(rho)+op_div(mx,my);
projerr(:,:,1) = projerr(:,:,1) - sizsc(1)*rho0;
projerr = norm(projerr(:),'inf');
% disp(norm(projerr(:),'inf'))
% disp(norm(projerr(:),2)/sqrt(prod(siz)))

%% 
    function divm = op_div(mx,my)
%         mx = m{1}; my = m{2};
        divm = sizsc(1)*cat(1,mx(1,:,:),...
                              mx(2:end,:,:)-mx(1:end-1,:,:),...
                             -mx(end,:,:))...
              +sizsc(2)*cat(2,my(:,1,:),...
                              my(:,2:end,:)-my(:,1:end-1,:),...
                             -my(:,end,:));
    end

    function [divadjphix,divadjphiy] = op_divadj(phi)
%         divadjphi = cell(1,2);
        divadjphix = sizsc(1)*(phi(1:end-1,:,:)-phi(2:end,:,:));
        divadjphiy = sizsc(2)*(phi(:,1:end-1,:)-phi(:,2:end,:));
    end

    function dtrho = op_dt(rho)
        dtrho = sizsc(3)*cat(3,rho(:,:,1),...
                               rho(:,:,2:end)-rho(:,:,1:end-1));
    end

    function dtadjphi = op_dtadj(phi)
        dtadjphi = sizsc(3)*cat(3,phi(:,:,1:end-1)-phi(:,:,2:end),...
                                  phi(:,:,end));
    end

end