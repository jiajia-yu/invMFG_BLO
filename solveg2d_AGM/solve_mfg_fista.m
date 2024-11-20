function [rho,mx,my,outs] = solve_mfg_fista(rho0,rho1,g,lambdaF,lambdaG,opts)
%% parameters
if isfield(opts,'nt') nt = opts.nt; else nt = 16; end
[nx,ny] = size(rho0);
siz = [nx,ny,nt];
dx = 1/nx; dy = 1/ny; dsp = dx*dy; dt = 1/nt;

% inital value of rho and m
if isfield(opts,'rho') xrho_nitm = opts.rho; else xrho_nitm = ones(nx,ny,nt); end
if isfield(opts,'mx')  xmx_nitm  = opts.mx;  else xmx_nitm  = zeros(nx-1,ny,nt); end
if isfield(opts,'my')  xmy_nitm  = opts.my;  else xmy_nitm  = zeros(nx,ny-1,nt); end

% stop criteria
if isfield(opts,'maxit') maxit = opts.maxit; else maxit = 5e3;  end
if isfield(opts,'tol')   tol   = opts.tol;   else tol   = 1e-6; end

% parameter for backtracking
if isfield(opts,'stepsize0') stepsize0 = opts.stepsize0; else stepsize0 = 0.1; end
if isfield(opts,'stepmodif') stepmodif = opts.stepmodif; else stepmodif = 0.8; end 
if isfield(opts,'submaxit')  submaxit  = opts.submaxit;  else submaxit  = 5; end

if isfield(opts,'acc') acc = logical(opts.acc); else acc = true(1); end

dotsc = @(a,b) sum(a.*b,'all')*dsp*dt;

%% initialization
obj_array = zeros(maxit,1);
res_array = zeros(maxit,1);
projerr_array = zeros(maxit,1);
stepsize_array = zeros(maxit,1);

yrho = xrho_nitm;
ymx = xmx_nitm;
ymy = xmy_nitm;
stepsize = stepsize0;
w_nitm = 1;
%% main iteration
for nit = 1:maxit
    % backtracking
    [objy,gradyrho,gradymx,gradymy] = comp_valgradL(siz,yrho,ymx,ymy,rho0,rho1,g,lambdaF,lambdaG);
    subnit = 1;
%     stepsize = stepsize0;
    while subnit < submaxit
        [xrho_nit,xmx_nit,xmy_nit,projerr] = proj_rhom(yrho - stepsize*gradyrho, ...
            ymx-stepsize*gradymx,ymy-stepsize*gradymy,rho0);
        xrho_nit = max(xrho_nit,1e-6);

        obj_nit = comp_valL(siz,xrho_nit,xmx_nit,xmy_nit,rho0,rho1,g,lambdaF,lambdaG);
        diffrho = xrho_nit - yrho;
        diffmx = xmx_nit - ymx;
        diffmy = xmy_nit - ymy;
        G = objy + dotsc(gradyrho,diffrho) + ...
                   dotsc(gradymx,diffmx)   + ...
                   dotsc(gradymy,diffmy)   + ...
            1/(2*stepsize)*(dotsc(diffrho,diffrho) + ...
                            dotsc(diffmx, diffmx)  + ...
                            dotsc(diffmy,diffmy));
        
        if obj_nit <= G
            break
        end
        subnit = subnit + 1;
        stepsize = stepsize*stepmodif;
    end
        
    % update variables
    w_nit = (1 + sqrt(1+4*w_nitm^2))/2;
    w = acc*(w_nitm-1)/w_nit;
    drho = xrho_nit-xrho_nitm;
    dmx = xmx_nit-xmx_nitm;
    dmy = xmy_nit-xmy_nitm;
    yrho = max( xrho_nit + w*drho,0.01);
%     yrho = xrhoNit + w*drho;
    ymx = xmx_nit + w*dmx;
    ymy = xmy_nit + w*dmy;
    
    w_nitm = w_nit;
    obj_nitm = obj_nit;
    xrho_nitm = xrho_nit;
    xmx_nitm = xmx_nit;
    xmy_nitm = xmy_nit;
    stepsize = max(stepsize,1e-5);
    
    res = sqrt(dotsc(drho,drho)+dotsc(dmx,dmx)+dotsc(dmy,dmy));
    
    obj_array(nit) = obj_nit;
    res_array(nit) = res;
    projerr_array(nit) = projerr;
    stepsize_array(nit) = stepsize;
    
    if res < tol
        break
    end
    fprintf('nit%d: %d sub iters, proj err %e, obj %f\n',...
             nit,subnit,projerr,obj_array(nit));
    
    if isnan(obj_nitm) || isinf(obj_nitm)
        fprintf('blow up at iteration %d with residue %.2e\n',nit,res);
        rho = xrho_nit;
        mx = xmx_nit;
        my = xmy_nit;
        outs.obj_array = obj_array(1:nit);
        outs.res_array = res_array(1:nit);
        outs.projerr_array = projerr_array(1:nit);
        outs.stepsize_array = stepsize_array(1:nit);
        return
    end
end

%% copy results
rho = xrho_nit;
mx = xmx_nit;
my = xmy_nit;
outs.obj_array = obj_array(1:nit);
outs.res_array = res_array(1:nit);
outs.projerr_array = projerr_array(1:nit);
outs.stepsize_array = stepsize_array(1:nit);

end