% main_test1
clear
clc
close all
rng(1)
% flg_gendata = true(1);
flg_gendata = false(1);

% metshape = 'qua';
metshape = 'wave';
nt = 16;
nx = 64; 

filename = [metshape,'_',num2str(nt),num2str(nx),];
if flg_gendata
    gendata(metshape,lambdab,nt,nx,nx);
else
    load(['../data/',filename,'.mat'])
end

%% solve BLO for g
g_min = min(g_true);
g_true = proj_g(g_true,g_min);
gtrue_norm = sqrt(mean(g_true.^2));
% parameters
T = 5e3;
K = 10;
% N = length(rhotilde_cell);
N = 1;
rhotol = 0.1;
valD_hist = zeros(T,1);
diffg_hist = zeros(T,1);
g_hist = zeros(size(g_true,1),T+1);

rho0_cell = rho0_cell(1:N);
rhotilde_cell = rhotilde_cell(1:N);
mtilde_cell = mtilde_cell(1:N,:);
siz = size(rhotilde_cell{1});
dd = 1/prod(siz);
siz1 = size(rhotilde_cell{1},1);
alpha = 8e0;
opts_initaux = [];
opts_initaux.betainit = 0.1;
opts_initaux.betaadj = 0.5;
opts_initaux.submaxit = 5;
opts_initaux.maxit = K;

lambdaL2 = 0;
mask_know = false(siz(1),1);
mask_know(1) = true;

% initialization
x = linspace(-1/2,1/2,siz1+1)';
x = (x(1:end-1)+x(2:end))/2;
g_num = mean(g_true)*ones(siz1,1);
g_num(mask_know) = g_true(mask_know);
mu_cell = rhotilde_cell;
w_cell = mtilde_cell;
g_hist(:,1) = g_num;

tic
% main iteration
for t = 1:T
    rho_cell = mu_cell;
    m_cell = w_cell;
    
    % initial auxiliary
    [valD_array,betak_array,rho_cellhist,m_cellhist] = initaux(mu_cell,w_cell,...
        rho0_cell,rho1_cell,g_num,lambdaF,lambdaG,rhotilde_cell,mtilde_cell,opts_initaux);
    
    valD_hist(t) = comp_valD(dd,rho_cellhist(:,end),m_cellhist(:,end),rhotilde_cell,mtilde_cell);
    valD_hist(t) = valD_hist(t) + lambdaL2*comp_valR_l2(g_num);
    
    % backtracking
    [gradg_bt,gradmu_cell,gradw_cell] = comp_backtrack(K,rho_cellhist,m_cellhist,...
        betak_array,mu_cell,w_cell,rho0_cell,rho1_cell,g_num,lambdaF,lambdaG,rhotilde_cell,mtilde_cell);
    gradg_l2 = comp_gradR_l2(g_num);
    gradg = gradg_bt + lambdaL2*gradg_l2;
    
    % update upper level variable
    g_num = proj_g(g_num - alpha*gradg,g_min);
    g_num(mask_know) = g_true(mask_know);
    diffg_hist(t) = sqrt(mean((g_num(:)-g_true(:)).^2))/gtrue_norm;
    g_hist(:,t+1) = g_num;
    
    mu_cell = rho_cellhist(:,end);
    w_cell = m_cellhist(:,end);
    fprintf('\nnum iter = %d, UL value = %f, rel err = %f\n',t,valD_hist(t),diffg_hist(t));
%     if mod(t,10)==0
%         tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
%         nexttile    
%         plot(x,g_true,'b','LineWidth',1.5);hold on;plot(x,g_num,'r','LineWidth',1.5);legend('true','num');
%         title(['T=',num2str(t)]);
%         pause(0.1)
%     end

end
t_run = toc;
filename = [metshape,'_',num2str(lambdaL2,'%.0e'),'_N=',num2str(N),'_T=',num2str(T),'_K=',num2str(K)];
save(['results/',filename]);

%%
fig=tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
nexttile    
plot(x,g_true,'b','LineWidth',1.5);hold on;plot(x,g_num,'r','LineWidth',1.5);legend('true','num');
%title(['T=',num2str(t)]);
exportgraphics(fig,['results/',filename,'_numg.png'],'BackgroundColor','none')

fig=tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
nexttile
plot(valD_hist);title('UL Objective');xlabel('UL iteration');ylabel('UL objective');
exportgraphics(fig,['results/',filename,'_ulobj.png'],'BackgroundColor','none')

fig=tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
nexttile
plot(diffg_hist,'linewidth',2);%title('||g_n-g_t||');xlabel('UL iteration');
exportgraphics(fig,['results/',filename,'_diffg.png'],'BackgroundColor','none')



