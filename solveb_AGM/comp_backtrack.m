function [gradb,gradmu_cell,gradw_cell] ...
    = comp_backtrack(kstar,rho_cellhist,m_cellhist,beta_array,...
    mu_cell,w_cell,rho0_cell,rho1_cell,b_num,lambdaF,lambdaG,rhotilde_cell,mtilde_cell)

N = length(rhotilde_cell);
siz = size(rhotilde_cell{1});
dt = 1/siz(3);
% dsp = 1/siz(1)/siz(2);
tol = 1e-8;
rhotol = 1e-1;

% initialization
gradb = zeros(siz(1,2));
gradmu_cell = cell(N,1);
gradw_cell = cell(N,2);

% last step
for n_data = 1:N
    gradmu_cell{n_data} = rho_cellhist{n_data,kstar} - rhotilde_cell{n_data};
    gradw_cell{n_data,1} = m_cellhist{n_data,1,kstar}- mtilde_cell{n_data,1};
    gradw_cell{n_data,2} = m_cellhist{n_data,2,kstar}- mtilde_cell{n_data,2};
end

% backtracking
if kstar == 1
    
    for n_data = 1:N
        gradb = gradb + op_Jb(gradmu_cell{n_data},...
                         {gradw_cell{n_data,1},gradw_cell{n_data,2}},...
                          beta_array(kstar));
        [gradmu,gradw] ...
            =op_Jmuw(rho0_cell{n_data},mu_cell{n_data},{w_cell{n_data,1},w_cell{n_data,2}},...
                 gradmu_cell{n_data},{gradw_cell{n_data,1},gradw_cell{n_data,2}},beta_array(1));
        gradmu_cell{n_data} = gradmu;
        gradw_cell{n_data,1} = gradw{1};
        gradw_cell{n_data,2} = gradw{2};
    end
    
else
    
    for k = kstar:-1:2
        for n_data = 1:N
        gradb = gradb + op_Jb(gradmu_cell{n_data},...
                         {gradw_cell{n_data,1},gradw_cell{n_data,2}},...
                          beta_array(kstar));
        [gradmu,gradw] ...
            =op_Jmuw(rho0_cell{n_data},rho_cellhist{n_data,k-1},{m_cellhist{n_data,1,k-1},m_cellhist{n_data,2,k-1}},...
                 gradmu_cell{n_data},{gradw_cell{n_data,1},gradw_cell{n_data,2}},beta_array(k));
        gradmu_cell{n_data} = gradmu;
        gradw_cell{n_data,1} = gradw{1};
        gradw_cell{n_data,2} = gradw{2};
            
        end
    end
    
    for n_data = 1:N
        gradb = gradb + op_Jb(gradmu_cell{n_data},...
                         {gradw_cell{n_data,1},gradw_cell{n_data,2}},...
                          beta_array(kstar));
        [gradmu,gradw] ...
            =op_Jmuw(rho0_cell{n_data},mu_cell{n_data},{w_cell{n_data,1},w_cell{n_data,2}},...
                 gradmu_cell{n_data},{gradw_cell{n_data,1},gradw_cell{n_data,2}},beta_array(1));
        gradmu_cell{n_data} = gradmu;
        gradw_cell{n_data,1} = gradw{1};
        gradw_cell{n_data,2} = gradw{2};        
    end
    
end


%% subfuncs

    function [gradmu_km,gradw_km] = op_Jmuw(rho0,rho,m,gradmu_k,gradw_k,beta)
        gradmu_temp = zeros(size(gradmu_k));
        gradwx_temp = zeros(size(gradw_k{1}));
        gradwy_temp = zeros(size(gradw_k{2}));

        gradmu_k(rho<rhotol) = 0;
        [gradmu_k,gradwx_k,gradwy_k] = proj_rhom(gradmu_k,gradw_k{1},gradw_k{2});
        gradmu_km = gradmu_k;
        gradw_k = {gradwx_k,gradwy_k};
        gradw_km = gradw_k;
        
        [Hrhorho,Hrhom,Hmm] = comp_dynHessian(rho0,rho,m);
        
        % gradmu
        gradmu_temp = gradmu_temp +0.25*dt*Hrhorho.*gradmu_k;
        gradmu_temp(:,:,1:end-1) = gradmu_temp(:,:,1:end-1) ...
            +dt*lambdaF./rho(:,:,1:end-1).*gradmu_k(:,:,1:end-1);
        gradmu_temp(:,:,1:end-1) = gradmu_temp(:,:,1:end-1) ...
            +0.25*dt*Hrhorho(:,:,2:end).*gradmu_k(:,:,1:end-1);
        gradmu_temp(:,:,end) = gradmu_temp(:,:,end) ...
            +lambdaG./rho(:,:,end).*gradmu_k(:,:,end);
        gradmu_temp(:,:,1:end-1) = gradmu_temp(:,:,1:end-1) ...
            +0.25*dt*Hrhorho(:,:,2:end).*gradmu_k(:,:,2:end);
        gradmu_temp(:,:,2:end) = gradmu_temp(:,:,2:end) ...
            +0.25*dt*Hrhorho(:,:,1:end-1).*gradmu_k(:,:,1:end-1);
        
        gradmu_temp(1:end-1,:,:) = gradmu_temp(1:end-1,:,:) ...
            +0.25*dt*Hrhom{1}(1:end-1,:,:).*gradw_k{1};
        gradmu_temp(2:end,:,:) = gradmu_temp(2:end,:,:) ...
            +0.25*dt*Hrhom{1}(2:end,:,:).*gradw_k{1};
        gradmu_temp(1:end-1,:,1:end-1) = gradmu_temp(1:end-1,:,1:end-1) ...
            +0.25*dt*Hrhom{1}(1:end-1,:,2:end).*gradw_k{1}(:,:,2:end);
        gradmu_temp(2:end,:,1:end-1) = gradmu_temp(2:end,:,1:end-1) ...
            +0.25*dt*Hrhom{1}(2:end,:,2:end).*gradw_k{1}(:,:,2:end);
        
        gradmu_temp(:,1:end-1,:) = gradmu_temp(:,1:end-1,:) ...
            +0.25*dt*Hrhom{2}(:,1:end-1,:).*gradw_k{2};
        gradmu_temp(:,2:end,:) = gradmu_temp(:,2:end,:) ...
            +0.25*dt*Hrhom{2}(:,2:end,:).*gradw_k{2};
        gradmu_temp(:,1:end-1,1:end-1) = gradmu_temp(:,1:end-1,1:end-1) ...
            +0.25*dt*Hrhom{2}(:,1:end-1,2:end).*gradw_k{2}(:,:,2:end);
        gradmu_temp(:,2:end,1:end-1) = gradmu_temp(:,2:end,1:end-1) ...
            +0.25*dt*Hrhom{2}(:,2:end,2:end).*gradw_k{2}(:,:,2:end);
        
        gradmu_km = gradmu_km - beta*gradmu_temp;
        
        % gradm
        gradwx_temp = gradwx_temp +0.25*dt*Hmm(1:end-1,:,:).*gradw_k{1};
        gradwx_temp = gradwx_temp +0.25*dt*Hmm(2:end,:,:).*gradw_k{1};
        gradwx_temp(1:end-1,:,:) = gradwx_temp(1:end-1,:,:) ...
            +0.25*dt*Hmm(2:end-1,:,:).*gradw_k{1}(2:end,:,:);
        gradwx_temp(2:end,:,:) = gradwx_temp(2:end,:,:) ...
            +0.25*dt*Hmm(2:end-1,:,:).*gradw_k{1}(1:end-1,:,:);
        
        gradwx_temp = gradwx_temp ...
            +0.25*dt*Hrhom{1}(1:end-1,:,:).*gradmu_k(1:end-1,:,:);
        gradwx_temp = gradwx_temp ...
            +0.25*dt*Hrhom{1}(2:end,:,:).*gradmu_k(2:end,:,:);
        gradwx_temp(:,:,2:end) = gradwx_temp(:,:,2:end) ...
            +0.25*dt*Hrhom{1}(1:end-1,:,2:end).*gradmu_k(1:end-1,:,1:end-1);
        gradwx_temp(:,:,2:end) = gradwx_temp(:,:,2:end) ...
            +0.25*dt*Hrhom{1}(2:end,:,2:end).*gradmu_k(2:end,:,1:end-1);
        
        gradwy_temp = gradwy_temp +0.25*dt*Hmm(:,1:end-1,:).*gradw_k{2};
        gradwy_temp = gradwy_temp +0.25*dt*Hmm(:,2:end,:).*gradw_k{2};
        gradwy_temp(:,1:end-1,:) = gradwy_temp(:,1:end-1,:) ...
            +0.25*dt*Hmm(:,2:end-1,:).*gradw_k{2}(:,2:end,:);
        gradwy_temp(:,2:end,:) = gradwy_temp(:,2:end,:) ...
            +0.25*dt*Hmm(:,2:end-1,:).*gradw_k{2}(:,1:end-1,:);
        
        gradwy_temp = gradwy_temp ...
            +0.25*dt*Hrhom{2}(:,1:end-1,:).*gradmu_k(:,1:end-1,:);
        gradwy_temp = gradwy_temp ...
            +0.25*dt*Hrhom{2}(:,2:end,:).*gradmu_k(:,2:end,:);
        gradwy_temp(:,:,2:end) = gradwy_temp(:,:,2:end) ...
            +0.25*dt*Hrhom{2}(:,1:end-1,2:end).*gradmu_k(:,1:end-1,1:end-1);
        gradwy_temp(:,:,2:end) = gradwy_temp(:,:,2:end) ...
            +0.25*dt*Hrhom{2}(:,2:end,2:end).*gradmu_k(:,2:end,1:end-1);
        
        gradw_km{1} = gradw_km{1} - beta*gradwx_temp;
        gradw_km{2} = gradw_km{2} - beta*gradwy_temp;
        
        
    end

    function term_km = op_Jb(gradmu_k,gradw_k,beta)
        [gradmu_k,~,~] = proj_rhom(gradmu_k,gradw_k{1},gradw_k{2});
        term_km = -beta*dt*sum(gradmu_k(:,:,1:end-1),3);
    end

    function [Hrhorho,Hrhom,Hmm] = comp_dynHessian(rho0,rho,m)
        mc = cell(1,2);
        mc{1} = cat(1, m{1}(1,:,:)/2, (m{1}(1:end-1,:,:)+m{1}(2:end,:,:))/2, m{1}(end,:,:)/2);
        mc{2} = cat(2, m{2}(:,1,:)/2, (m{2}(:,1:end-1,:)+m{2}(:,2:end,:))/2, m{2}(:,end,:)/2);
        msq = mc{1}.^2 + mc{2}.^2;
        rhoc = cat(3, (rho0+rho(:,:,1))/2, (rho(:,:,1:end-1)+rho(:,:,2:end))/2);
        
        ind = (rhoc>tol);
        Hrhorho = zeros(size(rhoc));
        Hrhorho(ind) = msq(ind)./rhoc(ind).^3;        
        Hmm = zeros(size(rhoc));
        Hmm(ind) = 1./rhoc(ind);
        Hrhom = cell(1,2);
        Hrhom{1} = zeros(size(rhoc));
        Hrhom{1}(ind) = -mc{1}(ind)./rhoc(ind).^2;
        Hrhom{2} = zeros(size(rhoc));
        Hrhom{2}(ind) = -mc{2}(ind)./rhoc(ind).^2;
        
    end

end