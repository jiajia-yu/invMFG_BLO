%% load or generate data

b_true = proj_b(b_true);
dd = 1/numel(rhotilde_cell{1});
% parameters
T = 1e3;
K = 10;
% N = length(rhotilde_cell);
N = 5;
rhotol = 0.1;
t_hist = zeros(T/10,1);
valD_hist = zeros(T/10,1);
diffb_hist = zeros(T,1);

rho0_cell = rho0_cell(1:N);
rhotilde_cell = rhotilde_cell(1:N);
mtilde_cell = mtilde_cell(1:N,:);
siz = size(rhotilde_cell{1});
alpha = 0.01;
% gamma = 0.1;
opts_initaux = [];
opts_initaux.betainit = 0.1;
opts_initaux.betaadj = 0.5;
opts_initaux.submaxit = 5;
opts_initaux.maxit = K;

% initialization
% b_num = zeros(siz(1:2));
b_num = rand(siz(1:2));
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
        rho0_cell,rho1,b_num,lambdaF,lambdaG,rhotilde_cell,mtilde_cell,opts_initaux);
    
%     valD_hist(t) = valD_array(end);
    
    % backtracking
    [gradb,gradmu_cell,gradw_cell] = comp_backtrack(K,rho_cellhist,m_cellhist,...
        betak_array,mu_cell,w_cell,rho0_cell,rho1,b_num,lambdaF,lambdaG,rhotilde_cell,mtilde_cell);
    
    % update upper level variable
    b_num = proj_b(b_num - alpha*gradb);
    diffb_hist(t) = sqrt(mean(((b_num(:)-b_true(:)).^2)));
    
    mu_cell = rho_cellhist(:,end);
    w_cell = m_cellhist(:,:,end);
    fprintf('\nnum iter = %d\n',t)
%     toc
    if mod(t,10)==0
%         t_hist(t/10) = t;
%         rho_tmp = cell(N,1);
%         m_tmp = cell(N,2);
%         for n_data = 1:N
%             [rho_tmp{n_data},m_tmp{n_data,1},m_tmp{n_data,2},outs] = ...
%                 solve_mfg_fista(rho0_cell{n_data},rho1,b_num,lambdaF,lambdaG,opts_gendata);
%         end
%         valD_hist(t/10) = comp_valD(dd,rho_tmp,m_tmp,rhotilde_cell,mtilde_cell);
        tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
        nexttile    
        imshow(b_num,[]);title(['numerical b, T=',num2str(t)]);colormap default; colorbar
        pause(0.1)
        
    end

end
t_run = toc;
save(['results/maze2r_N=',num2str(N),'_T=',num2str(t),'_K=',num2str(K)]);

%%
fig=tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
nexttile    
imshow(b_num,[]);title(['numerical b, T=',num2str(t),' K=',num2str(K)]);colormap default; colorbar
exportgraphics(fig,['results/maze2r_N=',num2str(N),'_T=',num2str(t),'_K=',num2str(K),'_numb.png'],'BackgroundColor','none')

fig=tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
nexttile
plot(t_hist(1:t/10),valD_hist(1:t/10));title('UL Objective (every 10 iters)');xlabel('UL iteration');ylabel('UL objective');
exportgraphics(fig,['results/maze2r_N=',num2str(N),'_T=',num2str(t),' K=',num2str(K),'_ulobj.png'],'BackgroundColor','none')

fig=tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
nexttile
plot(diffb_hist(1:t));title('||b_n-b_t||');xlabel('UL iteration');
exportgraphics(fig,['results/maze2r_N=',num2str(N),'_T=',num2str(t),' K=',num2str(K),'_diffb.png'],'BackgroundColor','none')

