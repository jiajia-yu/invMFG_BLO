% main_convergence
clear
clc
close all
% flg_gendata = true(1);
flg_gendata = false(1);

obsshape = 'gaussian';
lambdab = 0.05;
nt = 16;
nx = 64; 

filename = [obsshape,num2str(lambdab,'%.0e'),'_',num2str(nt),num2str(nx)];
if flg_gendata
    gendata(obsshape,lambdab,nt,nx,nx);
else
    load(['../data/',filename,'.mat'])
end

opts_gendata = [];
opts_gendata.maxit = 3e3;
opts_gendata.tol = 1e-6;
opts_gendata.stepsize0 = 0.01;
opts_gendata.stepmodif = 0.8;
opts_gendata.submaxit = 5;
opts_gendata.acc = true(1);

%% solve BLO for b
b_true = proj_b(b_true);
bnorm_true = sqrt(mean(b_true(:).^2));
dd = 1/numel(rhotilde_cell{1});
% parameters
T = 6e3;
K = 5;
% N = length(rhotilde_cell);
N = 1;
rhotol = 0.5;
t_hist = zeros(T/10+1,1);
b_hist = cell(T/10+1,1);
valD_hist = zeros(T/10+1,1);
diffb_rel_hist = zeros(T,1);

rho0_cell = rho0_cell(1:N);
rhotilde_cell = rhotilde_cell(1:N);
mtilde_cell = mtilde_cell(1:N,:);
siz = size(rhotilde_cell{1});
alpha = 1e-2;
gamma = 0.1;
opts_initaux = [];
opts_initaux.betainit = 0.1;
opts_initaux.betaadj = 0.5;
opts_initaux.submaxit = 5;
opts_initaux.maxit = K;

% initialization
b_num = zeros(siz(1:2));
% b_num = rand(siz(1:2));
b_best = b_num;
ind_best = 0;
diffb_rel_best = sqrt(mean((b_num(:)-b_true(:)).^2))/bnorm_true;
mu_cell = rhotilde_cell;
w_cell = mtilde_cell;

count_t = 1;
t_hist(count_t) = 0;
b_hist{count_t} = b_num;
valD_hist(count_t) = comp_valD(dd,mu_cell,w_cell,rhotilde_cell,mtilde_cell);

tic
% main iteration
for t = 1:T
    rho_cell = mu_cell;
    m_cell = w_cell;
    % initial auxiliary
    [valD_array,betak_array,rho_cellhist,m_cellhist] = initaux(mu_cell,w_cell,...
        rho0_cell,rho1_cell,b_num,lambdaF,lambdaG,rhotilde_cell,mtilde_cell,opts_initaux);
    
%     valD_hist(t) = valD_array(end);
    
    % backtracking
    [gradb,gradmu_cell,gradw_cell] = comp_backtrack(K,rho_cellhist,m_cellhist,...
        betak_array,mu_cell,w_cell,rho0_cell,rho1_cell,b_num,lambdaF,lambdaG,rhotilde_cell,mtilde_cell);
    
    % update upper level variable
    b_num = proj_b(b_num - alpha*gradb);
    diffb_rel = sqrt(mean((b_num(:)-b_true(:)).^2))/bnorm_true;
    diffb_rel_hist(t) = diffb_rel;
    if diffb_rel < diffb_rel_best
        diffb_rel_best = diffb_rel;
        b_best = b_num;
        ind_best = t;
    end
    
    mu_cell = rho_cellhist(:,end);
    w_cell = m_cellhist(:,:,end);
    fprintf('\nnum iter = %d\n',t)
%     toc
    if mod(t,10)==0
        count_t = count_t + 1;
        t_hist(count_t) = t;
        b_hist{count_t} = b_num;
        valD_hist(count_t) = comp_valD(dd,mu_cell,w_cell,rhotilde_cell,mtilde_cell);
%         tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
%         nexttile    
%         imshow(b_num,[]);title(['numerical b, T=',num2str(t)]);colormap default; colorbar
%         pause(0.1)
        
    end



end
t_run = toc;
% diffb_rel_hist = diffb_rel_hist./sqrt(mean(b_true(:).^2));
save(['results/',filename,'_new_N=',num2str(N),'_T=',num2str(T),'_K=',num2str(K)]);

%%
fig=tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
nexttile    
imshow(b_num,[]);title(['numerical b, T=',num2str(t),' K=',num2str(K)]);colormap default; colorbar
exportgraphics(fig,['results/',filename,'_N=',num2str(N),'_T=',num2str(T),'_K=',num2str(K),'_numb.png'],'BackgroundColor','none')

fig=tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
nexttile
plot(t_hist,valD_hist,'linewidth',2);
title('UL Objective');xlabel('UL iteration');ylabel('UL objective');
exportgraphics(fig,['results/',filename,'_N=',num2str(N),'_T=',num2str(T),'_K=',num2str(K),'_ulobj.png'],'BackgroundColor','none')

fig=tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
nexttile
plot(diffb_rel_hist,'linewidth',2);
title('obs relative error');xlabel('UL iteration');
exportgraphics(fig,['results/',filename,'_N=',num2str(N),'_T=',num2str(T),'_K=',num2str(K),'_diffb.png'],'BackgroundColor','none')



