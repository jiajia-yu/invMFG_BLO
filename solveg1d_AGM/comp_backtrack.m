function [gradg,gradmu_cell,gradw_cell] ...
    = comp_backtrack(kstar,rho_cellhist,m_cellhist,beta_array,...
    mu_cell,w_cell,rho0_cell,rho1_cell,g_num,lambdaF,lambdaG,rhotilde_cell,mtilde_cell)

N = length(rhotilde_cell);
siz = size(rhotilde_cell{1});
dt = 1/siz(2);
% dsp = 1/siz(1)/siz(2);
tol = 1e-8;
rhotol = 1e-1;
g_num = repmat(g_num,1,siz(2));

% initialization
gradg = zeros(siz(1),1);
gradmu_cell = cell(N,1);
gradw_cell = cell(N,1);

% last step
for n_data = 1:N
    gradmu_cell{n_data} = rho_cellhist{n_data,kstar} - rhotilde_cell{n_data};
    gradw_cell{n_data} = m_cellhist{n_data,kstar}- mtilde_cell{n_data};
end

% backtracking
if kstar == 1
    
    for n_data = 1:N
        gradg = gradg + op_Jb(rho0_cell{n_data},mu_cell{n_data},w_cell{n_data},...
                 gradmu_cell{n_data},gradw_cell{n_data,1},beta_array(kstar));
        [gradmu,gradw] ...
            =op_Jmuw(rho0_cell{n_data},mu_cell{n_data},w_cell{n_data},...
                 gradmu_cell{n_data},gradw_cell{n_data},beta_array(1));
        gradmu_cell{n_data} = gradmu;
        gradw_cell{n_data} = gradw;
    end
    
else
    
    for k = kstar:-1:2
        for n_data = 1:N
        gradg = gradg + op_Jb(rho0_cell{n_data},rho_cellhist{n_data,k-1},m_cellhist{n_data,k-1},...
                 gradmu_cell{n_data},gradw_cell{n_data},beta_array(k));
        [gradmu,gradw] ...
            =op_Jmuw(rho0_cell{n_data},rho_cellhist{n_data,k-1},m_cellhist{n_data,k-1},...
                 gradmu_cell{n_data},gradw_cell{n_data},beta_array(k));
        gradmu_cell{n_data} = gradmu;
        gradw_cell{n_data} = gradw;
            
        end
    end
    
    for n_data = 1:N
        gradg = gradg + op_Jb(rho0_cell{n_data},mu_cell{n_data},w_cell{n_data},...
                 gradmu_cell{n_data},gradw_cell{n_data,1},beta_array(1));
        [gradmu,gradw] ...
            =op_Jmuw(rho0_cell{n_data},mu_cell{n_data},w_cell{n_data},...
                 gradmu_cell{n_data},gradw_cell{n_data},beta_array(1));
        gradmu_cell{n_data} = gradmu;
        gradw_cell{n_data} = gradw;
    end
    
end


%% subfuncs

    function [gradmu_km,gradwx_km] = op_Jmuw(rho0,rho,m,gradmu_k,gradw_k,beta)
        gradmu_temp = zeros(size(gradmu_k));
        gradwx_temp = zeros(size(gradw_k));

        gradmu_k(rho<rhotol) = 0;
        [gradmu_k,gradwx_k] = proj_rhom(gradmu_k,gradw_k);
        gradmu_km = gradmu_k;
        gradw_k = gradwx_k;
        gradwx_km = gradw_k;
        
        [Hrhorho,Hrhom,Hmm] = comp_dynHessian(rho0,rho,m);
        
        % gradmu
        gradmu_temp = gradmu_temp +0.25*dt*Hrhorho.*gradmu_k;
        gradmu_temp(:,1:end-1) = gradmu_temp(:,1:end-1) ...
            +dt*lambdaF./rho(:,1:end-1).*gradmu_k(:,1:end-1);
        gradmu_temp(:,1:end-1) = gradmu_temp(:,1:end-1) ...
            +0.25*dt*Hrhorho(:,2:end).*gradmu_k(:,1:end-1);
        gradmu_temp(:,end) = gradmu_temp(:,end) ...
            +lambdaG./rho(:,end).*gradmu_k(:,end);
        gradmu_temp(:,1:end-1) = gradmu_temp(:,1:end-1) ...
            +0.25*dt*Hrhorho(:,2:end).*gradmu_k(:,2:end);
        gradmu_temp(:,2:end) = gradmu_temp(:,2:end) ...
            +0.25*dt*Hrhorho(:,1:end-1).*gradmu_k(:,1:end-1);
        
        gradmu_temp(1:end-1,:) = gradmu_temp(1:end-1,:) ...
            +0.25*dt*Hrhom(1:end-1,:).*gradw_k;
        gradmu_temp(2:end,:) = gradmu_temp(2:end,:) ...
            +0.25*dt*Hrhom(2:end,:).*gradw_k;
        gradmu_temp(1:end-1,1:end-1) = gradmu_temp(1:end-1,1:end-1) ...
            +0.25*dt*Hrhom(1:end-1,2:end).*gradw_k(:,2:end);
        gradmu_temp(2:end,1:end-1) = gradmu_temp(2:end,1:end-1) ...
            +0.25*dt*Hrhom(2:end,2:end).*gradw_k(:,2:end);
        
        gradmu_km = gradmu_km - beta*gradmu_temp;
        
        % gradm
        gradwx_temp = gradwx_temp +0.25*dt*Hmm(1:end-1,:).*gradw_k;
        gradwx_temp = gradwx_temp +0.25*dt*Hmm(2:end,:).*gradw_k;
        gradwx_temp(1:end-1,:) = gradwx_temp(1:end-1,:) ...
            +0.25*dt*Hmm(2:end-1,:).*gradw_k(2:end,:);
        gradwx_temp(2:end,:) = gradwx_temp(2:end,:) ...
            +0.25*dt*Hmm(2:end-1,:).*gradw_k(1:end-1,:);
        
        gradwx_temp = gradwx_temp ...
            +0.25*dt*Hrhom(1:end-1,:).*gradmu_k(1:end-1,:);
        gradwx_temp = gradwx_temp ...
            +0.25*dt*Hrhom(2:end,:).*gradmu_k(2:end,:);
        gradwx_temp(:,2:end) = gradwx_temp(:,2:end) ...
            +0.25*dt*Hrhom(1:end-1,2:end).*gradmu_k(1:end-1,1:end-1);
        gradwx_temp(:,2:end) = gradwx_temp(:,2:end) ...
            +0.25*dt*Hrhom(2:end,2:end).*gradmu_k(2:end,1:end-1);
                
        gradwx_km = gradwx_km - beta*gradwx_temp;
        
        
    end

    function term_km = op_Jb(rho0,rho_km,m_km,gradmu_k,gradw_k,beta)
        [gradmu_k,gradw_k,~] = proj_rhom(gradmu_k,gradw_k);
        [gradrho,gradm] = comp_dynGrad(rho0,rho_km,m_km);
        gradg_temp = 0.5*dt*sum(gradmu_k.*gradrho,2);
        gradg_temp(1:end-1) = gradg_temp(1:end-1) ...
                   + 0.5*dt*sum(gradmu_k(1:end-1,:).*gradrho(2:end,:),2);
        gradg_temp(1:end-1) = gradg_temp(1:end-1) ...
                   + 0.5*dt*sum(gradw_k.*gradm(1:end-1,:),2);
        gradg_temp(2:end) = gradg_temp(2:end) ...
                   + 0.5*dt*sum(gradw_k.*gradm(2:end,:),2);       
        term_km = -beta*gradg_temp;
    end

    function [gradrho,gradm] = comp_dynGrad(rho0,rho,m)
        mxc = cat(1, m(1,:)/2, (m(1:end-1,:)+m(2:end,:))/2, m(end,:)/2);
        msq = mxc.^2;
        rhoc = cat(2, (rho0+rho(:,1))/2, (rho(:,1:end-1)+rho(:,2:end))/2);
        
        ind = (rhoc>tol);
        gradrho = zeros(size(rhoc));
        gradrho(ind) = -msq(ind)./(2*rhoc(ind).^2);
        gradm = zeros(size(mxc));
        gradm(ind) = mxc(ind)./rhoc(ind);        
    end

    function [Hrhorho,Hrhomx,Hmm] = comp_dynHessian(rho0,rho,m)
        mxc = cat(1, m(1,:)/2, (m(1:end-1,:)+m(2:end,:))/2, m(end,:)/2);
        msq = mxc.^2;
        rhoc = cat(2, (rho0+rho(:,1))/2, (rho(:,1:end-1)+rho(:,2:end))/2);
        
        ind = (rhoc>tol);
        Hrhorho = zeros(size(rhoc));
        Hrhorho(ind) = g_num(ind).*msq(ind)./rhoc(ind).^3;        
        Hmm = zeros(size(rhoc));
        Hmm(ind) = g_num(ind)./rhoc(ind);
        Hrhomx = zeros(size(rhoc));
        Hrhomx(ind) = -g_num(ind).*mxc(ind)./rhoc(ind).^2;
        
    end

end