% test_hippocampus

%function 
% close all;
clear;

meshName = 'homer';
[surf.pt,surf.trg] = readOFF([meshName,'.off']);

nPt = size(surf.pt,1);
surf = surfOperators(surf);

%% example
% aroundpt = getaroundpt(surf.pt,surf.trg);
rho0 = zeros(nPt,1);
% rho0(surf.pt(:,2)<0.14 & surf.pt(:,2)>=0 & surf.pt(:,1)<0.1 & surf.pt(:,1)>-0.1) = 1;
% rho0(surf.pt(:,2)<0 & surf.pt(:,2)>-0.1) = 1;
rho0(surf.pt(:,2)>0.4) = 1;

rho1 = zeros(nPt,1);
% rho1(surf.pt(:,2)>0.4) = 1;
% rho1(surf.pt(:,3)>0.09 & surf.pt(:,2)>0 & surf.pt(:,2)<0.15) = 1;
rho1(surf.pt(:,2)<-0.43) = 1;
rho1 = rho1.*(sum(rho0.*surf.ptArea)/sum(rho1.*surf.ptArea));

rho0 = rho0 + 0.1;
rho1 = rho1 + 0.1;
logrho1 = log(rho1);

figure(1);clf(1);viewMesh(surf,rho0);colorbar
figure(2);clf(2);viewMesh(surf,rho1);colorbar

b = zeros(nPt,1);
b(surf.pt(:,3)>0 & surf.pt(:,2)>0 & surf.pt(:,2)<0.14...
                 & surf.pt(:,1)>0 & surf.pt(:,1)<0.12) = 1;
figure(3);clf(3);viewMesh(surf,b);axis on% colorbar

%% parameters
opts.funcL = @(rho,m) sum(m.^2,3)./(2*rho);
opts.gradLrho = @(rho,m) -sum(m.^2,3)./(2*rho.^2);
opts.gradLm = @(rho,m) m./rho; 
egName = 'obs';
lambdab = 1e1;
lambdae = 1e0;
opts.funcF = @(rho) lambdab*rho.*b + lambdae*rho.*log(rho);
opts.gradF = @(rho) lambdab*b + lambdae*(log(rho)+1);
lambdaG = 1e-1;
opts.funcG = @(rhoend) lambdaG*rhoend.*(log(rhoend)-logrho1);
opts.gradG = @(rhoend) lambdaG*(log(rhoend)-logrho1+1);

opts.plot = 1;
opts.savegif = 1;
opts.saveshot = 1;
opts.nt = 32;

opts.maxit = 3000;
opts.tol = 1e-7;

opts.stepsize0 = 4e-3;
opts.stepmodif = 0.5;
opts.submaxit = 5;
opts.acc = 1;


%% FISTA

tic;
% [rho, flux, output] = mfpMfFista(surf,rho0,rho1,opts);
[rho, flux, output] = mfgMfFista(surf,rho0,opts);
totaltime = toc
totalnit = length(output.objArray)
ittime = totaltime/totalnit

rhomin = min(rho,[],'all');
rhomax = max(rho,[],'all');

rhos = surf.pt2trg*(rho(:,1:end-1)+rho(:,2:end))/2;
% cost = output.costArray{end};
fprintf('dynamic cost: %f \n',...
        sum(surf.trgArea.*opts.funcL(rhos,flux),'all')/opts.nt);
fprintf('interaction cost: %f \n',...
        sum(surf.ptArea.*opts.funcF(rho(:,2:end-1)),'all')/opts.nt);
fprintf('terminal cost: %f \n',...
        sum(surf.ptArea.*opts.funcG(rho(:,end)),'all'));
fprintf('Total cost: %f \n',output.objArray(end));
figure;plot(sum(rho.*surf.ptArea));


filenameSave = [meshName,'_mfg_',egName];
save(['results/',filenameSave]);

%% visualization
% close all
if opts.plot
    %-------show and save objective history
    %h = figure(4);
    figure(1);
    plot(output.objArray(10:end),'LineWidth',2);
    print('-dpng',['results/' filenameSave,'_obj.png']);
    
    %-------show evolution
    figure(2);
%     ha = tight_subplot(1,2,[.01 .03],[.1 .01],[.01 .01]);
    set(gcf,'color','w');
%     set(gca,'position',[0 0 1 1],'units','normalized')
    for i = 1:size(rho,2)-1
        clf
        viewMesh(surf,rho(:,i));
        caxis([min(rho(:,i)),max(rho(:,i))]);colorbar;hold on
%         caxis([zmin,zmax]);colorbar;hold on
%         viewVectF(surf.trgCenter,squeeze(flux(:,i,:)));
        frame = getframe(gcf);
        im{i} = frame2im(frame);
        pause(0.2);        
    end
    i = size(rho,2);
    clf
    viewMesh(surf,rho(:,i));
    caxis([min(rho(:,i)),max(rho(:,i))]);colorbar;hold on
%     caxis([zmin,zmax]);colorbar;hold on
    frame = getframe(gcf);
    im{i} = frame2im(frame);
    pause(0.2);
    
    %--------- save gif
    if opts.savegif
        for idx = 1:size(rho,2)
            [A,map] = rgb2ind(im{idx},256);
            if idx == 1
                imwrite(A,map,['results/',filenameSave,'.gif'],'gif','LoopCount',Inf,'DelayTime',0.3);
            else
                imwrite(A,map,['results/',filenameSave,'.gif'],'gif','WriteMode','append','DelayTime',0.3);
            end
        end
    end

    %---------- save eps
    if opts.saveshot
        if ~isfield(opts,'num_frame') num_frame = 5; end
        idx_frame = round(linspace(1,opts.nt+1,num_frame));
        t_frame = (idx_frame-1)./(opts.nt);

        close all
        for k = 1:num_frame-1
            clf
            idx = idx_frame(k);
            viewMesh(surf,rho(:,idx));hold on;
%             viewVectF(surf.trgCenter(1:5:end,:),squeeze(flux(1:5:end,idx,:)));
            set(gcf,'unit','centimeters','position',[10 5 2 3])
            set(gca,'Position',[0.1,0.05,0.65,0.8]);% left margin, lower margin, width, height
            caxis([min(rho(:,idx)),max(rho(:,idx))]);
%             caxis([zmin,zmax]);colorbar;hold on
            colorbar('Position',[0.76,0.1,0.05,0.7]);% left margin, lower margin, width,
            title(['t=',num2str(t_frame(k))]);
            fig = gcf;
            exportgraphics(fig,['results/',filenameSave,'_shot',num2str(k),'.eps']);
        end
        clf
        k = num_frame;
        idx = idx_frame(k);
        viewMesh(surf,rho(:,idx));hold on;
%         viewVectF(surf.trgCenter(1:5:end,:),squeeze(flux(1:5:end,idx,:)));
        set(gcf,'unit','centimeters','position',[10 5 2 3])
        set(gca,'Position',[0.1,0.05,0.65,0.8]);% left margin, lower margin, width, height
        caxis([min(rho(:,idx)),max(rho(:,idx))]);
%         caxis([zmin,zmax]);colorbar;hold on
        colorbar('Position',[0.76,0.1,0.05,0.7]);% left margin, lower margin, width,
        title(['t=',num2str(t_frame(k))]);
        fig = gcf;
        exportgraphics(fig,['results/',filenameSave,'_shot',num2str(k),'.eps']);
        
    end
    
end

