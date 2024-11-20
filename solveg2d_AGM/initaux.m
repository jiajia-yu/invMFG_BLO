function [valD_array,beta_array,rho_cellhist,m_cellhist] = initaux(mu_cell,w_cell,...
    rho0_cell,rho1_cell,g,lambdaF,lambdaG,rhotilde_cell,mtilde_cell,opts)

if isfield(opts,'betainit') beta_init = opts.betainit; else, beta_init = 0.1; end
if isfield(opts,'betaadj') beta_adj = opts.betaadj; else, beta_adj = 0.5; end
if isfield(opts,'submaxit') submaxit = opts.submaxit; else, submaxit = 5; end
if isfield(opts,'maxit') maxit = opts.maxit; else, maxit = 20; end
rhotol = 1e-1;

N = length(mu_cell);
siz = size(mu_cell{1});
dd = 1/prod(siz);
beta_array = zeros(N,maxit);
valD_array = zeros(maxit,1);
obj_array = zeros(N,1);
rho_cellhist = cell(N,maxit);
m_cellhist = cell(N,2,maxit);
dotsc = @(a,b) sum(a.*b,'all')*dd;

rho_cell = mu_cell;
m_cell = w_cell;

for n_data = 1:N
    rho0 = rho0_cell{n_data};
    rho1 = rho1_cell{n_data};
    rho_nit = rho_cell{n_data};
    mx_nit = m_cell{n_data,1};
    my_nit = m_cell{n_data,2};
    obj_array(n_data) = comp_valL(siz,rho_nit,mx_nit,my_nit,rho0,rho1,g,lambdaF,lambdaG);
end

% main iteration
for nit = 1:maxit
    for n_data = 1:N
        rho0 = rho0_cell{n_data};
        rho1 = rho1_cell{n_data};
        rho_nit = rho_cell{n_data};
        mx_nit = m_cell{n_data,1};
        my_nit = m_cell{n_data,2};
        obj_nit = obj_array(n_data);
        % backtracking
        [gradyrho,gradymx,gradymy] = comp_gradL(siz,rho_nit,mx_nit,my_nit,rho0,rho1,g,lambdaF,lambdaG);
        subnit = 1;
        beta_nit = beta_init;
        while subnit < submaxit
            [rho_nitp,mx_nitp,my_nitp,projerr] = proj_rhom(rho_nit - beta_nit*gradyrho, ...
                mx_nit-beta_nit*gradymx,my_nit-beta_nit*gradymy,rho0);
            rho_nitp = max(rho_nitp,rhotol);

            obj_nitp = comp_valL(siz,rho_nitp,mx_nitp,my_nitp,rho0,rho1,g,lambdaF,lambdaG);
            diffrho = rho_nitp - rho_nit;
            diffmx = mx_nitp - mx_nit;
            diffmy = my_nitp - my_nit;
            G = obj_nit + dotsc(gradyrho,diffrho) + ...
                       dotsc(gradymx,diffmx)   + ...
                       dotsc(gradymy,diffmy)   + ...
                1/(2*beta_nit)*(dotsc(diffrho,diffrho) + ...
                                dotsc(diffmx,diffmx)  + ...
                                dotsc(diffmy,diffmy));

            if obj_nitp <= G
                break
            end
            subnit = subnit + 1;
            beta_nit = beta_nit*beta_adj;
        end

        % update variables
        obj_array(n_data) = obj_nitp;
        rho_cell{n_data} = rho_nitp;
        rho_cellhist{n_data,nit} = rho_nitp;
        m_cell{n_data,1} = mx_nitp;
        m_cell{n_data,2} = my_nitp;
        m_cellhist{n_data,1,nit} = mx_nitp;
        m_cellhist{n_data,2,nit} = my_nitp;
        beta_array(n_data,nit) = beta_nit;

%         fprintf('nit%d: %d sub iters, proj err %e, obj %f\n',...
%                  nit,subnit,projerr,obj_array(n_data));
%              
%         if isnan(obj_nit) || isinf(obj_nit)
%             fprintf('blow up at iteration %d with residue %.2e\n',nit,res);
%             rho = cat(2,rho0_cell,xrhoNit);
%             m = xmNit;
%             outs.objArray = objArray(1:nit);
%             outs.resArray = resArray(1:nit);
%             outs.costArray = costArray(1:nit);
%             outs.projerrArray = projerrArray(1:nit);
%             outs.stepsizeArray = stepsizeArray(1:nit);
%             return
%         end
    end
    valD_array(nit) = comp_valD(dd,rho_cell,m_cell,rhotilde_cell,mtilde_cell);
    
end

end