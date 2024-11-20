load('results\2dtri_1e-04_N=4_T=5000_K=5.mat')

fprintf('diffg=%f, time=%f\n',diffg_hist(end),t_run);

%%
fig=tiledlayout(1,3,'TileSpacing','Compact','Padding','Compact');
set(gcf,'unit','centimeters','Position',[9 9 16 5]);
nexttile    
imshow(g_true+4,[]);colormap default; colorbar
nexttile    
imshow(g_true+2,[]);colormap default; colorbar
nexttile    
imshow(g_true+1,[]);colormap default; colorbar
exportgraphics(fig,['results/',filename,'_trueg.eps'],'BackgroundColor','none')

fig=tiledlayout(1,3,'TileSpacing','Compact','Padding','Compact');
set(gcf,'unit','centimeters','Position',[9 9 16 5]);
nexttile    
imshow(g_num{1},[]);colormap default; colorbar
nexttile    
imshow(g_num{2},[]);colormap default; colorbar
nexttile    
imshow(g_num{3},[]);colormap default; colorbar
exportgraphics(fig,['results/',filename,'_numg.eps'],'BackgroundColor','none')

fig=tiledlayout(1,3,'TileSpacing','Compact','Padding','Compact');
set(gcf,'unit','centimeters','Position',[9 9 16 5]);
nexttile    
imshow(g_num{1}-(g_true+4),[]);colormap default; colorbar
nexttile    
imshow(g_num{2}-(g_true+2),[]);colormap default; colorbar
nexttile    
imshow(g_num{3}-(g_true+1),[]);colormap default; colorbar
exportgraphics(fig,['results/',filename,'_errg.eps'],'BackgroundColor','none')

%%
fig=tiledlayout(1,3,'TileSpacing','Compact','Padding','Compact');
set(gcf,'unit','centimeters','Position',[9 9 16 5]);
nexttile    
imshow(g_true+4,[]);colormap default; colorbar
nexttile    
imshow(g_num{1},[]);colormap default; colorbar
nexttile    
imshow(g_num{1}-(g_true+4),[]);colormap default; colorbar
exportgraphics(fig,['results/',filename,'_gxx.eps'],'BackgroundColor','none')

fig=tiledlayout(1,3,'TileSpacing','Compact','Padding','Compact');
set(gcf,'unit','centimeters','Position',[9 9 16 5]);
nexttile    
imshow(g_true+2,[]);colormap default; colorbar
nexttile    
imshow(g_num{2},[]);colormap default; colorbar
nexttile    
imshow(g_num{2}-(g_true+2),[]);colormap default; colorbar
exportgraphics(fig,['results/',filename,'_gxy.eps'],'BackgroundColor','none')

fig=tiledlayout(1,3,'TileSpacing','Compact','Padding','Compact');
set(gcf,'unit','centimeters','Position',[9 9 16 5]);
nexttile    
imshow(g_true+1,[]);colormap default; colorbar
nexttile    
imshow(g_num{3},[]);colormap default; colorbar
nexttile    
imshow(g_num{3}-(g_true+1),[]);colormap default; colorbar
exportgraphics(fig,['results/',filename,'_gyy.eps'],'BackgroundColor','none')

