%% convergence - data
clear
% filename = 'gaussian1e+00_1664';
% filename = 'gaussian5e-01_1664';
% filename = 'gaussian1e-01_1664';
filename = 'gaussian5e-02_1664';
% filename = 'ring1e-01_3264';
load(['../data/',filename,'.mat'])
ind = [1,4,8,12,16];

for ind_i = 1:5
    i = ind(ind_i);
    rho = squeeze(rhotilde_cell{1}(:,:,i));
    mx = squeeze(mtilde_cell{1,1}(:,:,i));
    mx = cat(1,mx(1,:)/2,(mx(1:end-1,:)+mx(2:end,:))/2,mx(end,:)/2);
    my = squeeze(mtilde_cell{1,2}(:,:,i));
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
    exportgraphics(fig,['rev/',filename,'_data',num2str(i),'.eps'],'BackgroundColor','none')
end

%% robustness - result
clear
filename = 'twobar5e-01_uniform0e+00_N=1_T=5000_K=5';
% filename = 'twobar5e-01_uniform3e-01_N=1_T=5000_K=5';
% filename = 'twobar5e-01_uniform5e-01_N=1_T=5000_K=5';
% filename = 'twobar5e-01_uniform8e-01_N=1_T=5000_K=5';
load(['results/',filename,'.mat'])
ind = [1,4,8,12,16];

for ind_i = 1:5
    i = ind(ind_i);
    rho = squeeze(rhotilde_cell{1}(:,:,i));
    mx = squeeze(mtilde_cell{1,1}(:,:,i));
    mx = cat(1,mx(1,:)/2,(mx(1:end-1,:)+mx(2:end,:))/2,mx(end,:)/2);
    my = squeeze(mtilde_cell{1,2}(:,:,i));
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
    exportgraphics(fig,['rev/',filename,'_data',num2str(i),'.eps'],'BackgroundColor','none')
end

%% varobs - result
clear
% filename = 'ring5e-01_N=1_T=5000_K=5';
filename = 'clover5e-01_N=1_T=5000_K=5';
load(['results/',filename,'.mat'])
ind = [1,8,16,24,32];

for ind_i = 1:5
    i = ind(ind_i);
    rho = squeeze(rhotilde_cell{1}(:,:,i));
    mx = squeeze(mtilde_cell{1,1}(:,:,i));
    mx = cat(1,mx(1,:)/2,(mx(1:end-1,:)+mx(2:end,:))/2,mx(end,:)/2);
    my = squeeze(mtilde_cell{1,2}(:,:,i));
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
    exportgraphics(fig,['rev/',filename,'_data',num2str(i),'.eps'],'BackgroundColor','none')
end
