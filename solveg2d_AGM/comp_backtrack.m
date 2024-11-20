function [gradg,gradmu_cell,gradw_cell] ...
    = comp_backtrack(kstar,rho_cellhist,m_cellhist,beta_array,...
    mu_cell,w_cell,rho0_cell,rho1_cell,g_num,lambdaF,lambdaG,rhotilde_cell,mtilde_cell)

N = length(rhotilde_cell);
siz = size(rhotilde_cell{1});
dt = 1/siz(3);
% dsp = 1/siz(1)/siz(2);
tol = 1e-8;
rhotol = 1e-1;
gxx = g_num{1};
gxy = g_num{2};
gyy = g_num{3};

% initialization
gradg = cell(3,1);
gradg{1} = zeros(siz(1,2));
gradg{2} = zeros(siz(1,2));
gradg{3} = zeros(siz(1,2));
gradmu_cell = cell(N,1);
gradw_cell = cell(N,2);

% last step
for n_data = 1:N
    gradmu_cell{n_data} = rho_cellhist{n_data,kstar} - rhotilde_cell{n_data};
    gradw_cell{n_data,1} = m_cellhist{n_data,1,kstar}- mtilde_cell{n_data,1};
    gradw_cell{n_data,2} = m_cellhist{n_data,2,kstar}- mtilde_cell{n_data,2};
end

% backtracking
    
for k = kstar:-1:2
    for n_data = 1:N
        [rhoc,mc] = stg2ctr(rho0_cell{n_data},rho_cellhist{n_data,k-1},{m_cellhist{n_data,1,k-1},m_cellhist{n_data,2,k-1}});
        gradmu_cell{n_data}(rho_cellhist{n_data,k-1}<rhotol) = 0;
        [gradmu_cell{n_data},gradw_cell{n_data,1},gradw_cell{n_data,2}] ...
            = proj_rhom(gradmu_cell{n_data},gradw_cell{n_data,1},gradw_cell{n_data,2});

        term_km = op_Jb(rhoc,mc,gradmu_cell{n_data},{gradw_cell{n_data,1},gradw_cell{n_data,2}},beta_array(k));
        gradg{1} = gradg{1} + term_km{1};
        gradg{2} = gradg{2} + term_km{2};
        gradg{3} = gradg{3} + term_km{3};
        [gradmu,gradw] = op_Jmuw(rho_cellhist{n_data,k-1},rhoc,mc,gradmu_cell{n_data},{gradw_cell{n_data,1},gradw_cell{n_data,2}},beta_array(k));
        gradmu_cell{n_data} = gradmu;
        gradw_cell{n_data,1} = gradw{1};
        gradw_cell{n_data,2} = gradw{2};

    end
end

for n_data = 1:N
    [rhoc,mc] = stg2ctr(rho0_cell{n_data},mu_cell{n_data},{w_cell{n_data,1},w_cell{n_data,2}});
    gradmu_cell{n_data}(mu_cell{n_data}<rhotol) = 0;
    [gradmu_cell{n_data},gradw_cell{n_data,1},gradw_cell{n_data,2}] ...
        = proj_rhom(gradmu_cell{n_data},gradw_cell{n_data,1},gradw_cell{n_data,2});
    
    term_km = op_Jb(rhoc,mc,gradmu_cell{n_data},{gradw_cell{n_data,1},gradw_cell{n_data,2}},beta_array(k));
    gradg{1} = gradg{1} + term_km{1};
    gradg{2} = gradg{2} + term_km{2};
    gradg{3} = gradg{3} + term_km{3};    
    [gradmu,gradw] = op_Jmuw(mu_cell{n_data},rhoc,mc,gradmu_cell{n_data},{gradw_cell{n_data,1},gradw_cell{n_data,2}},beta_array(1));
    gradmu_cell{n_data} = gradmu;
    gradw_cell{n_data,1} = gradw{1};
    gradw_cell{n_data,2} = gradw{2};        
end
    


%% subfuncs

    function [gradmu_km,gradw_km] = op_Jmuw(rho_km,rhoc,mc,gradmu_k,gradw_k,beta)
        gradmu_temp = zeros(size(gradmu_k));
        gradwx_temp = zeros(size(gradw_k{1}));
        gradwy_temp = zeros(size(gradw_k{2}));

        gradmu_km = gradmu_k;
        gradw_km = gradw_k;
        
        [Hrhorho,Hrhom,Hmm] = comp_dynHessian(rhoc,mc);
        gradmu_cent = 0.25*dt*cat(3, gradmu_k(:,:,1), gradmu_k(:,:,1:end-1)+gradmu_k(:,:,2:end));
        gradw_cent{1} = 0.25*dt*cat(1, gradw_k{1}(1,:,:), gradw_k{1}(1:end-1,:,:)+gradw_k{1}(2:end,:,:), gradw_k{1}(end,:,:));
        gradw_cent{2} = 0.25*dt*cat(2, gradw_k{2}(:,1,:), gradw_k{2}(:,1:end-1,:)+gradw_k{2}(:,2:end,:), gradw_k{2}(:,end,:));
        
        % gradmu
        % HrhorhoF, HrhorhoG
        gradmu_temp(:,:,1:end-1) = gradmu_temp(:,:,1:end-1) ...
            +dt*lambdaF./rho_km(:,:,1:end-1).*gradmu_k(:,:,1:end-1);
        gradmu_temp(:,:,end) = gradmu_temp(:,:,end) ...
            +lambdaG./rho_km(:,:,end).*gradmu_k(:,:,end);
        % HrhorhoL
        gradmu_temp = gradmu_temp +Hrhorho.*gradmu_cent;
        gradmu_temp(:,:,1:end-1) = gradmu_temp(:,:,1:end-1) ...
            +Hrhorho(:,:,2:end).*gradmu_cent(:,:,2:end);
        % HrhomxL
        gradmu_temp = gradmu_temp +Hrhom{1}.*gradw_cent{1};
        gradmu_temp(:,:,1:end-1) = gradmu_temp(:,:,1:end-1) ...
            +Hrhom{1}(:,:,2:end).*gradw_cent{1}(:,:,2:end);
        % HrhomyL
        gradmu_temp = gradmu_temp +Hrhom{2}.*gradw_cent{2};
        gradmu_temp(:,:,1:end-1) = gradmu_temp(:,:,1:end-1) ...
            +Hrhom{2}(:,:,2:end).*gradw_cent{2}(:,:,2:end);
        
        gradmu_km = gradmu_km - beta*gradmu_temp;
        
        % gradmx
        % Hrhomx
        gradwx_temp = gradwx_temp +Hrhom{1}(1:end-1,:,:).*gradmu_cent(1:end-1,:,:);
        gradwx_temp = gradwx_temp +Hrhom{1}(2:end,:,:).*gradmu_cent(2:end,:,:);
        % Hmxmx
        gradwx_temp = gradwx_temp +Hmm{1}(1:end-1,:,:).*gradw_cent{1}(1:end-1,:,:);
        gradwx_temp = gradwx_temp +Hmm{1}(2:end,:,:).*gradw_cent{1}(2:end,:,:);
       % Hmxmy
        gradwx_temp = gradwx_temp +Hmm{2}(1:end-1,:,:).*gradw_cent{2}(1:end-1,:,:);
        gradwx_temp = gradwx_temp +Hmm{2}(2:end,:,:).*gradw_cent{2}(2:end,:,:);
        % gradmy
        % Hrhomy
        gradwy_temp = gradwy_temp +Hrhom{2}(:,1:end-1,:).*gradmu_cent(:,1:end-1,:);
        gradwy_temp = gradwy_temp +Hrhom{2}(:,2:end,:).*gradmu_cent(:,2:end,:);
        % Hmxmy
        gradwy_temp = gradwy_temp +Hmm{2}(:,1:end-1,:).*gradw_cent{1}(:,1:end-1,:);
        gradwy_temp = gradwy_temp +Hmm{2}(:,2:end,:).*gradw_cent{1}(:,2:end,:);
        % Hmymy
        gradwy_temp = gradwy_temp +Hmm{3}(:,1:end-1,:).*gradw_cent{2}(:,1:end-1,:);
        gradwy_temp = gradwy_temp +Hmm{3}(:,2:end,:).*gradw_cent{2}(:,2:end,:);
        
        gradw_km{1} = gradw_km{1} - beta*gradwx_temp;
        gradw_km{2} = gradw_km{2} - beta*gradwy_temp;
        
        
    end

    function term_km = op_Jb(rhoc,mc,gradmu_k,gradw_k,beta)
        term_km = cell(3,1);
%         gradmu_k(rho_km<rhotol) = 0;
%         [gradmu_k,gradwx_k,gradwy_k] = proj_rhom(gradmu_k,gradw_k{1},gradw_k{2});
        gradmu_cent = 0.25*dt*cat(3, gradmu_k(:,:,1), gradmu_k(:,:,1:end-1)+gradmu_k(:,:,2:end));
        gradw_cent{1} = 0.25*dt*cat(1, gradw_k{1}(1,:,:), gradw_k{1}(1:end-1,:,:)+gradw_k{1}(2:end,:,:), gradw_k{1}(end,:,:));
        gradw_cent{2} = 0.25*dt*cat(2, gradw_k{2}(:,1,:), gradw_k{2}(:,1:end-1,:)+gradw_k{2}(:,2:end,:), gradw_k{2}(:,end,:));
        
        ind = (rhoc>tol);
        rhoinv = 1./rhoc;
        rhoinv(~ind) = 0;
        
        temp = sum(-(mc{1}.^2.*rhoinv.^2/2).*gradmu_cent,3);
        temp = temp + sum(mc{1}.*rhoinv.*gradw_cent{1},3);
        term_km{1} = -beta*dt*temp;
        temp = sum(-(mc{1}.*mc{2}.*rhoinv.^2/2).*gradmu_cent,3);
        temp = temp + sum(mc{2}.*rhoinv.*gradw_cent{1},3);
        temp = temp + sum(mc{1}.*rhoinv.*gradw_cent{2},3);
        term_km{2} = -beta*dt*temp;
        temp = sum(-(mc{2}.^2.*rhoinv.^2/2).*gradmu_cent,3);
        temp = temp + sum(mc{2}.*rhoinv.*gradw_cent{2},3);
        term_km{3} = -beta*dt*temp;
    end

    function [Hrhorho,Hrhom,Hmm] = comp_dynHessian(rhoc,mc)
%         mc = cell(1,2);
%         mc{1} = cat(1, m{1}(1,:,:)/2, (m{1}(1:end-1,:,:)+m{1}(2:end,:,:))/2, m{1}(end,:,:)/2);
%         mc{2} = cat(2, m{2}(:,1,:)/2, (m{2}(:,1:end-1,:)+m{2}(:,2:end,:))/2, m{2}(:,end,:)/2);
        msq = gxx.*mc{1}.^2 + gxy.*mc{1}.*mc{2} + gyy.*mc{2}.^2;
%         rhoc = cat(3, (rho0+rho(:,:,1))/2, (rho(:,:,1:end-1)+rho(:,:,2:end))/2);
        
        ind = (rhoc>tol);
        Hrhorho = zeros(size(rhoc));
        Hrhorho(ind) = msq(ind)./rhoc(ind).^3;        
        Hmm = cell(1,3);
        Hmm{1} = gxx./rhoc; Hmm{1}(~ind) = 0;
        Hmm{2} = gxy./rhoc; Hmm{2}(~ind) = 0;
        Hmm{3} = gyy./rhoc; Hmm{3}(~ind) = 0;
        Hrhom = cell(1,2);
        Hrhom{1} = -(gxx.*mc{1} + gxy.*mc{2})./rhoc.^2; Hrhom{1}(~ind) = 0;
        Hrhom{2} = -(gxy.*mc{1} + gyy.*mc{2})./rhoc.^2; Hrhom{2}(~ind) = 0;
        
    end

    function [rhoc,mc] = stg2ctr(rho0,rho,m)
        rhoc = cat(3, (rho0+rho(:,:,1))/2, (rho(:,:,1:end-1)+rho(:,:,2:end))/2);
        mc = cell(1,2);
        mc{1} = cat(1, m{1}(1,:,:)/2, (m{1}(1:end-1,:,:)+m{1}(2:end,:,:))/2, m{1}(end,:,:)/2);
        mc{2} = cat(2, m{2}(:,1,:)/2, (m{2}(:,1:end-1,:)+m{2}(:,2:end,:))/2, m{2}(:,end,:)/2);
        
    end

end