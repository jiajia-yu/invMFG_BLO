%% convergence - data
clear
% filename = 'gaussian1e+00_1664';
% filename = 'gaussian5e-01_1664';
% filename = 'gaussian1e-01_1664';
filename = 'gaussian5e-02_1664';
% filename = 'ring1e-01_3264';
load(['../data/',filename,'.mat'])

rho = squeeze(rhotilde_cell{1}(:,:,8));
mx = squeeze(mtilde_cell{1,1}(:,:,8));
mx = cat(1,mx(1,:)/2,(mx(1:end-1,:)+mx(2:end,:))/2,mx(end,:)/2);
my = squeeze(mtilde_cell{1,2}(:,:,8));
my = cat(2,my(:,1)/2,(my(:,1:end-1)+my(:,2:end))/2,my(:,end)/2);
% mc = sqrt(mx(32,32)^2+my(32,32)^2)
mnorm = sqrt(mx.^2 + my.^2);
mx(mnorm<prctile(mnorm(:),70)) = 0;
my(mnorm<prctile(mnorm(:),70)) = 0;
rho(32,32)

nx = 64; ny = 64;
vec_dens = 4;
vec_leng = 0.5;
X = repmat( (1:vec_dens:nx) ,floor(ny/vec_dens),1 );
Y = repmat( (1:vec_dens:ny)',1, floor(nx/vec_dens));

fig=tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
set(gcf,'unit','centimeters','Position',[9 9 7 5]);
nexttile    
imshow(rho,[]);colormap default; colorbar
colorbar('Position',[0.85,0.15,0.03,0.8]);
hold on;
quiver(X,Y,...
       my(1:vec_dens:end,1:vec_dens:end),mx(1:vec_dens:end,1:vec_dens:end),...
       vec_leng,'r','LineWidth',0.8);
exportgraphics(fig,['illu/',filename,'_data.eps'],'BackgroundColor','none')

%% convergence - result random init
clear
% filename = 'gaussian1e+00_1664_N=1_T=6000_K=5';
% filename = 'gaussian5e-01_1664_N=1_T=6000_K=5';
% filename = 'gaussian1e-01_1664_N=1_T=6000_K=5';
filename = 'gaussian5e-02_1664_N=1_T=6000_K=5';
%5e-02,1e-01,,5e+00
% filename = 'ring5e+00_3264_N=1_T=6000_K=5';
load(['results/',filename,'.mat'])

fprintf('stepsize = %f, ul value = %f\n', alpha, valD_hist(end));
fprintf('best ind = %d, diffb best = %f, diffb end=%f\n',...
         ind_best,diffb_rel_best,diffb_rel);

fig=tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
set(gcf,'unit','centimeters','Position',[9 9 7 5]);
nexttile    
imshow(b_best,[]);%title(['numerical b, T=',num2str(t),' K=',num2str(K)]);
colormap default; colorbar
exportgraphics(fig,['figresults/',filename,'_N=',num2str(N),'_T=',num2str(T),'_K=',num2str(K),'_numb.eps'],'BackgroundColor','none')

fig=tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
set(gcf,'unit','centimeters','Position',[9 9 7 5]);
nexttile
plot(t_hist,valD_hist,'linewidth',2);title('UL Objective');xlabel('UL iteration');ylabel('UL objective');
exportgraphics(fig,['figresults/',filename,'_N=',num2str(N),'_T=',num2str(T),'_K=',num2str(K),'_ulobj.eps'],'BackgroundColor','none')

fig=tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
set(gcf,'unit','centimeters','Position',[9 9 7 5]);
nexttile
plot(diffb_rel_hist,'linewidth',2);title('obs relative error');xlabel('UL iteration');
exportgraphics(fig,['figresults/',filename,'_N=',num2str(N),'_T=',num2str(T),'_K=',num2str(K),'_diffb.eps'],'BackgroundColor','none')

%% convergence - result zero init
clear
% filename = 'gaussian1e+00_1664_new_N=1_T=6000_K=5.mat';
% filename = 'gaussian5e-01_1664_new_N=1_T=6000_K=5.mat';
% filename = 'gaussian1e-01_1664_new_N=1_T=6000_K=5.mat';
filename = 'gaussian5e-02_1664_new_N=1_T=6000_K=5.mat';

load(['results\',filename])
% filename = 'gaussian1e+00_1664';
% filename = 'gaussian5e-01_1664';
% filename = 'gaussian1e-01_1664';
filename = 'gaussian5e-02_1664';

fprintf('stepsize = %f, ul value = %f\n', alpha, valD_hist(end));
fprintf('best ind = %d, diffb best = %f, diffb end=%f, time=%f\n',...
         ind_best,diffb_rel_best,diffb_rel,t_run);

fig=tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
set(gcf,'unit','centimeters','Position',[9 9 7 5]);
nexttile    
imshow(b_best,[]);%title(['numerical b, T=',num2str(t),' K=',num2str(K)]);
colormap default; colorbar
colorbar('Position',[0.85,0.15,0.03,0.8]);
exportgraphics(fig,['figresults/',filename,'_N=',num2str(N),'_T=',num2str(T),'_K=',num2str(K),'_numb.eps'],'BackgroundColor','none')

fig=tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
set(gcf,'unit','centimeters','Position',[9 9 7 5]);
nexttile
plot(t_hist,valD_hist,'linewidth',2);title('UL Objective');xlabel('UL iteration');ylabel('UL objective');
exportgraphics(fig,['figresults/',filename,'_N=',num2str(N),'_T=',num2str(T),'_K=',num2str(K),'_ulobj.eps'],'BackgroundColor','none')

fig=tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
set(gcf,'unit','centimeters','Position',[9 9 7 5]);
nexttile
plot(diffb_rel_hist,'linewidth',2);title('obs relative error');xlabel('UL iteration');
exportgraphics(fig,['figresults/',filename,'_N=',num2str(N),'_T=',num2str(T),'_K=',num2str(K),'_diffb.eps'],'BackgroundColor','none')

%% robustness - result
clear
filename = 'twobar5e-01_uniform0e+00_N=1_T=5000_K=5';
% filename = 'twobar5e-01_uniform3e-01_N=1_T=5000_K=5';
% filename = 'twobar5e-01_uniform5e-01_N=1_T=5000_K=5';
% filename = 'twobar5e-01_uniform8e-01_N=1_T=5000_K=5';
load(['results/',filename,'.mat'])

fprintf('best ind = %d, diffb best = %f, diffb end=%f, time=%f\n',...
         ind_best,diffb_rel_best,diffb_rel,t_run);

fig=tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
set(gcf,'unit','centimeters','Position',[9 9 7 5]);
nexttile    
imshow(squeeze(rhotilde_cell{1}(:,:,8)),[]);%title(['numerical b, T=',num2str(t),' K=',num2str(K)]);
colormap default; colorbar
colorbar('Position',[0.85,0.15,0.03,0.8]);
exportgraphics(fig,['figresults/',filename,'_N=',num2str(N),'_T=',num2str(T),'_K=',num2str(K),'_rho.eps'],'BackgroundColor','none')

fig=tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
set(gcf,'unit','centimeters','Position',[9 9 7 5]);
nexttile    
imshow(b_num,[]);%title(['numerical b, T=',num2str(t),' K=',num2str(K)]);
colormap default; colorbar
colorbar('Position',[0.85,0.15,0.03,0.8]);
exportgraphics(fig,['figresults/',filename,'_N=',num2str(N),'_T=',num2str(T),'_K=',num2str(K),'_numb.eps'],'BackgroundColor','none')

fig=tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
set(gcf,'unit','centimeters','Position',[9 9 7 5]);
nexttile    
imshow(b_num-b_true,[]);%title(['numerical b, T=',num2str(t),' K=',num2str(K)]);
colormap default; colorbar
colorbar('Position',[0.85,0.15,0.03,0.7]);
exportgraphics(fig,['figresults/',filename,'_N=',num2str(N),'_T=',num2str(T),'_K=',num2str(K),'_diffb.eps'],'BackgroundColor','none')

fig=tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
set(gcf,'unit','centimeters','Position',[9 9 7 5]);
nexttile
plot(diffb_rel_hist,'linewidth',2);title('obs relative error');xlabel('UL iteration');
exportgraphics(fig,['figresults/',filename,'_N=',num2str(N),'_T=',num2str(T),'_K=',num2str(K),'_diffbhist.eps'],'BackgroundColor','none')

%% varobs - result
clear
filename = 'ring5e-01_N=1_T=5000_K=5';
% filename = 'clover5e-01_N=1_T=5000_K=5';
load(['results/',filename,'.mat'])

fprintf('best ind = %d, diffb best = %f, diffb end=%f\n',...
         ind_best,diffb_rel_best,diffb_rel);

fig=tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
set(gcf,'unit','centimeters','Position',[9 9 7 5]);
nexttile    
imshow(b_true,[]);%title(['numerical b, T=',num2str(t),' K=',num2str(K)]);
colormap default; colorbar
colorbar('Position',[0.85,0.15,0.03,0.8]);
exportgraphics(fig,['figresults/',filename,'_N=',num2str(N),'_T=',num2str(T),'_K=',num2str(K),'_trueb.eps'],'BackgroundColor','none')

fig=tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
set(gcf,'unit','centimeters','Position',[9 9 7 5]);
nexttile    
imshow(b_best,[]);%title(['numerical b, T=',num2str(t),' K=',num2str(K)]);
colormap default; colorbar
colorbar('Position',[0.85,0.15,0.03,0.8]);
exportgraphics(fig,['figresults/',filename,'_N=',num2str(N),'_T=',num2str(T),'_K=',num2str(K),'_numb.eps'],'BackgroundColor','none')

fig=tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
set(gcf,'unit','centimeters','Position',[9 9 7 5]);
nexttile    
imshow(b_best-b_true,[]);%title(['numerical b, T=',num2str(t),' K=',num2str(K)]);
colormap default; colorbar
colorbar('Position',[0.85,0.15,0.03,0.8]);
exportgraphics(fig,['figresults/',filename,'_N=',num2str(N),'_T=',num2str(T),'_K=',num2str(K),'_diffb.eps'],'BackgroundColor','none')

fig=tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
set(gcf,'unit','centimeters','Position',[9 9 7 5]);
nexttile
plot(diffb_rel_hist,'linewidth',2);title('obs relative error');xlabel('UL iteration');
exportgraphics(fig,['figresults/',filename,'_N=',num2str(N),'_T=',num2str(T),'_K=',num2str(K),'_diffbhist.eps'],'BackgroundColor','none')
