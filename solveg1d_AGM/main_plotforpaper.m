%% robustness - result
clear
load('results\cub_uniform0e+001e-05_N=1_T=5000_K=10.mat')
% load('results\cub_uniform1e-013e-04_N=1_T=5000_K=10.mat')
% load('results\cub_uniform2e-011e-03_N=1_T=5000_K=10.mat')
% load('results\cub_uniform3e-013e-03_N=1_T=5000_K=10.mat')

fprintf('sigma = %f, lambda L2 = %f, diffg end=%f, time=%f\n',...
         sigma,lambdaL2,diffg_hist(end),t_run);

fig=tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
set(gcf,'unit','centimeters','Position',[9 9 7 5]);
nexttile    
[T,X] = meshgrid(linspace(0,1,nt+1),x);
mesh(T,X,cat(2,rho0_cell{1},rhotilde_cell{1}));
xlabel('t');ylabel('x');
exportgraphics(fig,['figresults/',filename,'_datarho.eps'],'BackgroundColor','none')

fig=tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
set(gcf,'unit','centimeters','Position',[9 9 7 5]);
nexttile    
[T,X] = meshgrid(linspace(1/nt,1-1/nt,nt),(x(1:end-1)+x(2:end))/2);
mesh(T,X,mtilde_cell{1});
xlabel('t');ylabel('x');
exportgraphics(fig,['figresults/',filename,'_datam.eps'],'BackgroundColor','none')

fig=tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
set(gcf,'unit','centimeters','Position',[9 9 7 5]);
nexttile    
plot(x,g_true,'b','LineWidth',1.5);hold on;
plot(x,g_num,'r','LineWidth',1.5);legend('true','num','location','north');
ylim([0.5,1.5])
exportgraphics(fig,['figresults/',filename,'_numg.eps'],'BackgroundColor','none')

fig=tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
set(gcf,'unit','centimeters','Position',[9 9 7 5]);
nexttile    
plot(diffg_hist,'linewidth',2);%title('||g_n-g_t||');xlabel('UL iteration');
exportgraphics(fig,['figresults/',filename,'_diffg.eps'],'BackgroundColor','none')

%% varmets
load('results\qua_1e-05_N=1_T=5000_K=10.mat')

fprintf('lambda L2 = %f, diffg end=%f, time=%f\n',...
         lambdaL2,diffg_hist(end),t_run);

fig=tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
set(gcf,'unit','centimeters','Position',[9 9 7 5]);
nexttile    
plot(x,g_true,'b','LineWidth',1.5);hold on;
plot(x,g_num,'r','LineWidth',1.5);legend('true','num','location','south');
exportgraphics(fig,['figresults/',filename,'_numg.eps'],'BackgroundColor','none')

fig=tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
set(gcf,'unit','centimeters','Position',[9 9 7 5]);
nexttile    
plot(diffg_hist,'linewidth',2);%title('||g_n-g_t||');xlabel('UL iteration');
exportgraphics(fig,['figresults/',filename,'_diffg.eps'],'BackgroundColor','none')

%---
load('results\wave_1e-05_N=1_T=5000_K=10.mat')

fprintf('lambda L2 = %f, diffg end=%f, time=%f\n',...
         lambdaL2,diffg_hist(end),t_run);

fig=tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
set(gcf,'unit','centimeters','Position',[9 9 7 5]);
nexttile    
plot(x,g_true,'b','LineWidth',1.5);hold on;
plot(x,g_num,'r','LineWidth',1.5);
ylim([0.7,2.6]),legend('true','num','location','northeast');
exportgraphics(fig,['figresults/',filename,'_numg.eps'],'BackgroundColor','none')

fig=tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
set(gcf,'unit','centimeters','Position',[9 9 7 5]);
nexttile    
plot(diffg_hist,'linewidth',2);%title('||g_n-g_t||');xlabel('UL iteration');
exportgraphics(fig,['figresults/',filename,'_diffg.eps'],'BackgroundColor','none')

%% multidata
clear
load('results\tri_0e+00_N=1_T=5000_K=10.mat')
% load('results\tri_0e+00_N=2_T=5000_K=10.mat')
% load('results\tri_1e-05_N=1_T=5000_K=10.mat')
% load('results\tri_1e-04_N=2_T=5000_K=10.mat')

fprintf('lambda L2 = %f, diffg end=%f, time=%f\n',...
         lambdaL2,diffg_hist(end),t_run);

fig=tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
set(gcf,'unit','centimeters','Position',[9 9 7 5]);
nexttile    
plot(x,g_true,'b','LineWidth',1.5);hold on;
plot(x,g_num,'r','LineWidth',1.5);legend('true','num','location','north');
exportgraphics(fig,['figresults/',filename,'_numg.eps'],'BackgroundColor','none')

fig=tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
set(gcf,'unit','centimeters','Position',[9 9 7 5]);
nexttile    
plot(diffg_hist,'linewidth',2);%title('||g_n-g_t||');xlabel('UL iteration');
exportgraphics(fig,['figresults/',filename,'_diffg.eps'],'BackgroundColor','none')

