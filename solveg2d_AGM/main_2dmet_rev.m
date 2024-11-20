% main_test1
clear
clc
close all
% rng(1)
% flg_gendata = true(1);
flg_gendata = false(1);

metshape = '2dtri';
nt = 30;
nx = 50; 

filename = [metshape,'_',num2str(nt),num2str(nx),];
if flg_gendata
    gendata(nt,nx,nx);
else
    load(['../data/',filename,'.mat'])
end

%% solve BLO for g 
dd = 1/numel(rhotilde_cell{1});
gtol = comp_gtol({g_true+4,g_true+2,g_true+1});
gtrue_norm = sqrt(mean((g_true(:)+4).^2) ...
                + mean((g_true(:)+2).^2) ...
                + mean((2*g_true(:)+1).^2));
% opts_gendata = [];
% opts_gendata.maxit = 5e3;
% opts_gendata.tol = 1e-8;
% opts_gendata.stepsize0 = 0.001;
% opts_gendata.stepmodif = 0.8;
% opts_gendata.submaxit = 5;
% opts_gendata.acc = true(1);
% parameters
T = 5e3;
K = 5;
% N = length(rhotilde_cell);
N = 4;
rhotol = 0.1;
rho0_cell = rho0_cell(1:N);
rho1_cell = rho1_cell(1:N);
rhotilde_cell = rhotilde_cell(1:N);
mtilde_cell = mtilde_cell(1:N,:);
siz = size(rhotilde_cell{1});
step_uv = 1e2;
step_lv = 0.1;
% t_hist = zeros(T/10,1);
valD_hist = zeros(T,1);
diffg_hist = zeros(T,1);
lambdaL2 = 0;
lambdaP = 0;
mask = false(siz(1,2));
mask(1,:) = true;

opts_initaux = [];
opts_initaux.betainit = 0.1;
opts_initaux.betaadj = 0.5;
opts_initaux.submaxit = 5;
opts_initaux.maxit = K;

% initialization
g_num = {4*ones(siz(1:2));2*ones(siz(1:2));ones(siz(1:2))};
% g_num = {rand(siz(1:2));rand(siz(1:2));rand(siz(1:2))};
mu_cell = rhotilde_cell;
w_cell = mtilde_cell;
% mu_cell = cell(N,1);
% w_cell = cell(N,2);
% for n = 1:N
%     mu_cell{n} = ones(siz);
%     w_cell{n,1} = zeros(siz-[1,0,0]);
%     w_cell{n,2} = zeros(siz-[0,1,0]);
% end

tic
% main iteration
for t = 1:T
    rho_cell = mu_cell;
    m_cell = w_cell;
    
    % initial auxiliary
    [valD_array,betak_array,rho_cellhist,m_cellhist] = initaux(mu_cell,w_cell,...
        rho0_cell,rho1_cell,g_num,lambdaF,lambdaG,rhotilde_cell,mtilde_cell,opts_initaux);
    
    valD_hist(t) = valD_array(end);
    valD_hist(t) = valD_hist(t) + lambdaL2*comp_valR_l2(g_num);
    valD_hist(t) = valD_hist(t) + lambdaP*comp_valR_p(g_num,{g_true+4;g_true+2;g_true+1},mask);
    
    % backtracking
    [gradg_bt,gradmu_cell,gradw_cell] = comp_backtrack(K,rho_cellhist,m_cellhist,...
        betak_array,mu_cell,w_cell,rho0_cell,rho1_cell,g_num,lambdaF,lambdaG,rhotilde_cell,mtilde_cell);
    gradg_l2 = comp_gradR_l2(g_num);
    gradg_p = comp_gradR_p(g_num,{g_true+4;g_true+2;g_true+1},mask);
    
    % update upper level variable
    g_num{1} = g_num{1} - step_uv*(gradg_bt{1}+lambdaL2*gradg_l2{1}+lambdaP*gradg_p{1});
    g_num{2} = g_num{2} - step_uv*(gradg_bt{2}+lambdaL2*gradg_l2{2}+lambdaP*gradg_p{2});
    g_num{3} = g_num{3} - step_uv*(gradg_bt{3}+lambdaL2*gradg_l2{3}+lambdaP*gradg_p{3});
    g_num = proj_g(g_num,gtol);
    diffg_hist(t) = sqrt(mean((g_num{1}(:)-g_true(:)-4).^2) ...
                        +mean((g_num{2}(:)-g_true(:)-2).^2) ...
                        +mean((g_num{3}(:)-g_true(:)-1).^2))/gtrue_norm;
    
    mu_cell = rho_cellhist(:,end);
    w_cell = m_cellhist(:,:,end);
    fprintf('\nnum iter = %d, diffg = %e\n',t,diffg_hist(t))
%     toc
%     if mod(t,10)==0
% %         t_hist(t/10) = t;
% %         rho_tmp = cell(N,1);
% %         m_tmp = cell(N,2);
% %         for n_data = 1:N
% %             [rho_tmp{n_data},m_tmp{n_data,1},m_tmp{n_data,2},outs] = ...
% %                 solve_mfg_fista(rho0_cell{n_data},rho1,g_num,lambdaF,lambdaG,opts_gendata);
% %         end
% %         valD_hist(t/10) = comp_valD(dd,rho_tmp,m_tmp,rhotilde_cell,mtilde_cell);
%         clf;
%         plot(diffg_hist(1:t),'linewidth',2)
% 
% %         tiledlayout(2,2,'TileSpacing','Compact','Padding','Compact');
% %         nexttile    
% %         imshow(g_num{1},[]);title(['numerical gxx, T=',num2str(t)]);colormap default; colorbar
% %         nexttile    
% %         imshow(g_num{2},[]);title(['numerical gxy, T=',num2str(t)]);colormap default; colorbar
% %         nexttile    
% %         imshow(g_num{3},[]);title(['numerical gyy, T=',num2str(t)]);colormap default; colorbar
% %         nexttile    
% %         imshow(g_true,[]);title(['ground truth g, T=',num2str(t)]);colormap default; colorbar
%         pause(0.1)
%         
%     end

end
t_run = toc;
filename = [metshape,'_rev_',num2str(lambdaL2,'%.0e'),...
                        '_N=',num2str(N),'_T=',num2str(T),'_K=',num2str(K)];
save(['results/',filename]);

%%

fig=tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
nexttile    
imshow(g_true,[]);colormap default; colorbar
% title(['numerical gxx, T=',num2str(t)]);
% exportgraphics(fig,['results/',filename,'_trueg.png'],'BackgroundColor','none')

fig=tiledlayout(1,3,'TileSpacing','Compact','Padding','Compact');
nexttile    
imshow(g_num{1}-(g_true+4),[]);colormap default; colorbar
% title(['numerical gxx, T=',num2str(t)]);
nexttile    
imshow(g_num{2}-(g_true+2),[]);colormap default; colorbar
% title(['numerical gxy, T=',num2str(t)]);
nexttile    
imshow(g_num{3}-(g_true+1),[]);colormap default; colorbar
% title(['numerical gyy, T=',num2str(t)]);
% nexttile    
% imshow(g_true,[]);title(['ground truth g, T=',num2str(t)]);colormap default; colorbar
% exportgraphics(fig,['results/',filename,'_errg.png'],'BackgroundColor','none')

fig=tiledlayout(1,3,'TileSpacing','Compact','Padding','Compact');
nexttile    
imshow(g_num{1},[]);colormap default; colorbar
% title(['numerical gxx, T=',num2str(t)]);
nexttile    
imshow(g_num{2},[]);colormap default; colorbar
% title(['numerical gxy, T=',num2str(t)]);
nexttile    
imshow(g_num{3},[]);colormap default; colorbar
% title(['numerical gyy, T=',num2str(t)]);
% nexttile    
% imshow(g_true,[]);title(['ground truth g, T=',num2str(t)]);colormap default; colorbar
% exportgraphics(fig,['results/',filename,'_numg.png'],'BackgroundColor','none')


% fig=tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
% nexttile
% plot(valD_hist);title('UL Objective');xlabel('UL iteration');ylabel('UL objective');
% exportgraphics(fig,['results/',filename,'_ulobj.png'],'BackgroundColor','none')
% 
fig=tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
nexttile
plot(diffg_hist);title('||g_n-g_t||');xlabel('UL iteration');
% exportgraphics(fig,['results/',filename,'_diffg.png'],'BackgroundColor','none')



