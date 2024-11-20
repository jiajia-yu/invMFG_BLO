% main_test1
clear
clc
close all
% gendata = true(1);
gendata = false(1);

%% generate/load training data
if ~gendata
    load('../data/grav1d2.mat');
else
    N =1;
    % mesh
    nt = 16; ntp = nt+1; ntm = nt-1;
    dt = 1/nt;

    nx = 64; nxp = nx+1; nxm = nx-1;
    dx = 1/nx;

    x = reshape( linspace(-1/2,1/2,nxp),nxp,1);
    x = (x(1:end-1,:)+x(2:end,:))/2;

    % block, terminal and parameters
    eg = 'grav1d2';
    lambdaF = 1e-2;
    lambdaG = 5e-1;
%     rho1 = 10*exp(-( (x).^2 )*100 ) + 0.1;
    g_true = -0.3*cos(2*pi*x) + 0.7;
%     g_true = -cos(2*pi*x) + 2;
    
    % opts
    opts_gendata = [];
    opts_gendata.maxit = 5e3;
    opts_gendata.tol = 1e-6;
    opts_gendata.stepsize0 = 0.001;
    opts_gendata.stepmodif = 0.8;
    opts_gendata.submaxit = 5;
    opts_gendata.acc = true(1);

    % initial and computation
    initcen = [-0.3; -0.2; -0.1; 0; 0.1];
    rho0_cell = cell(N,1);
%     rho0_cell{1} = 0.8+0.4*(sin(2*pi*(x+0.5)).^2);
    rho0_cell{1} = 10*exp(-( (x).^2 )*100 ) + 0.1;
%     rho0_cell{2} = 10*exp(-((abs(x)-0.5).^2)*100)+0.1;
%     rho0_cell{3} = 10*exp(-( (x-0.3).^2 )*100 ) + 0.1;
%     rho0_cell{4} = 10*exp(-( (x+0.3).^2 )*100 ) + 0.1;
    
    rhotilde_cell = cell(N,1); 
    mtilde_cell = cell(N,1);
    for n_data = 1:N
        rho0 = rho0_cell{n_data};
%         rho1 = 0.8+0.4*(sin(pi*(x+0.5)).^2);
        rho1 = mean(rho0)*ones(size(x));
        
        tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
        nexttile
        plot(x,rho0,'b');hold on; plot(x,rho1,'r'); plot(x,g_true,'k');
        rho0_cell{n_data} = rho0;

        [rho,mx,outs] = solve_mfg_fista(rho0,rho1,g_true,lambdaF,lambdaG,opts_gendata);
        figure;plot(outs.obj_array);
        show_movement(cat(2,rho0,rho),rho1,[],['illu/grav1d2_',num2str(n_data)])
        rhotilde_cell{n_data} = rho;
        mtilde_cell{n_data} = mx;
    end
%     save('../data/grav1d2.mat','rho1','g_true','lambdaF','lambdaG','rho0_cell','rhotilde_cell','mtilde_cell');
end

%% perturb data
noisetype = 'uniform';
sigma = 0;
tol = 0;
switch noisetype
    case 'gaussian'
        for n = 1:size(rho0_cell,1)
            rho0_cell{n} = max(tol, rho0_cell{n}...
                            +sigma*randn(size(rho0_cell{n})));
            rhotilde_cell{n} = max(tol, rhotilde_cell{n}...
                            +sigma*randn(size(rhotilde_cell{n})));
            mtilde_cell{n,1} = mtilde_cell{n,1}...
                            +sigma*randn(size(mtilde_cell{n,1}));
        end
    case 'uniform'
        for n = 1:size(rho0_cell,1)
            rho0_cell{n} = max(tol, rho0_cell{n}...
                        +sigma*(rand(size(rho0_cell{n}))-0.5));
            rhotilde_cell{n} = max(tol, rhotilde_cell{n}...
                        +sigma*(rand(size(rhotilde_cell{n}))-0.5));
            mtilde_cell{n,1} = mtilde_cell{n,1}...
                        +sigma*(rand(size(mtilde_cell{n,1}))-0.5);
        end
end


%% solve BLO for b
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
alpha = 1e0;
opts_initaux = [];
opts_initaux.betainit = 0.1;
opts_initaux.betaadj = 0.5;
opts_initaux.submaxit = 5;
opts_initaux.maxit = K;

lambdaL2 = 1e-4;
lambdaP = 1e-2;
mask = false(siz(1),1);
% mask([1 siz1/2,end]) = true;
mask(1) = true;

% initialization
x = linspace(-1/2,1/2,siz1+1)';
x = (x(1:end-1)+x(2:end))/2;
% g_num = mean(g_true)*ones(siz(1),1);
g_num = 0.5*(max(g_true)+min(g_true))+ 0.5*(max(g_true)-min(g_true))*randn(siz1(1),1);
mu_cell = rhotilde_cell;
w_cell = mtilde_cell;
g_hist(:,1) = g_num;

% fig = tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
% nexttile    
% plot(x,g_true,'b','LineWidth',1.5);hold on;plot(x,g_num,'r','LineWidth',1.5);legend('true','num');
% title(['T=',num2str(0)]);
% pause(0.1)
% exportgraphics(fig,['videos/grav1d2_N=',num2str(N),'_K=',num2str(K),'_0.png'],'BackgroundColor','none')
% 
% count = 1;

tic
% main iteration
for t = 1:T
    rho_cell = mu_cell;
    m_cell = w_cell;
    
    % initial auxiliary
    [valD_array,betak_array,rho_cellhist,m_cellhist] = initaux(mu_cell,w_cell,...
        rho0_cell,rho1,g_num,lambdaF,lambdaG,rhotilde_cell,mtilde_cell,opts_initaux);
    
    valD_hist(t) = comp_valD(dd,rho_cellhist(:,end),m_cellhist(:,end),rhotilde_cell,mtilde_cell);
    valD_hist(t) = valD_hist(t) + lambdaL2*comp_valR_l2(g_num);
    valD_hist(t) = valD_hist(t) + lambdaP*comp_valR_p(g_num,g_true,mask);
    
    % backtracking
    [gradg_bt,gradmu_cell,gradw_cell] = comp_backtrack(K,rho_cellhist,m_cellhist,...
        betak_array,mu_cell,w_cell,rho0_cell,rho1,g_num,lambdaF,lambdaG,rhotilde_cell,mtilde_cell);
    gradg_l2 = comp_gradR_l2(g_num);
    gradg_p = comp_gradR_p(g_num,g_true,mask);
    gradg = gradg_bt + lambdaL2*gradg_l2 + lambdaP*gradg_p;
    
    % update upper level variable
    g_num = proj_g(g_num - alpha*gradg,g_min);
    diffg_hist(t) = sqrt(mean((g_num(:)-g_true(:)).^2))/gtrue_norm;
    g_hist(:,t+1) = g_num;
    
    mu_cell = rho_cellhist(:,end);
    w_cell = m_cellhist(:,end);
    fprintf('\nnum iter = %d, UL value = %f, rel err = %f\n',t,valD_hist(t),diffg_hist(t));
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
        plot(x,g_true,'b','LineWidth',1.5);hold on;plot(x,g_num,'r','LineWidth',1.5);legend('true','num');
        title(['T=',num2str(t)]);
        pause(0.1)
%         exportgraphics(fig,['videos/grav1d2_N=',num2str(N),'_K=',num2str(K),'_',num2str(count),'.png'],'BackgroundColor','none')
%         count = count + 1;
    end

end
t_run = toc;
filename = ['grav1d2_','_',noisetype,num2str(sigma,'%.0e'),num2str(lambdaL2,'%.0e'),'_',num2str(lambdaP,'%.0e'),...
            '_N=',num2str(N),'_T=',num2str(T),'_K=',num2str(K)];
% save(['results/',filename]);

%%
fig=tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
nexttile    
plot(x,g_true,'b','LineWidth',1.5);hold on;plot(x,g_num,'r','LineWidth',1.5);legend('true','num');
%title(['T=',num2str(t)]);
% exportgraphics(fig,['results/',filename,'_numg.png'],'BackgroundColor','none')

fig=tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
nexttile
plot(valD_hist);title('UL Objective');xlabel('UL iteration');ylabel('UL objective');
% exportgraphics(fig,['results/',filename,'_ulobj.png'],'BackgroundColor','none')

fig=tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
nexttile
plot(diffg_hist,'linewidth',2);%title('||g_n-g_t||');xlabel('UL iteration');
% exportgraphics(fig,['results/',filename,'_diffg.png'],'BackgroundColor','none')



