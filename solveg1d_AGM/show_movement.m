function show_movement(rho,rho1,opts,filename)
[nx,ntp] = size(rho);
x = reshape( linspace(-1/2,1/2,nx+1),nx+1,1);
x = (x(1:end-1,:)+x(2:end,:))/2;

% zmin = min(rho,[],'all');
% zmax = max(rho,[],'all');

if nargin < 3
    opts = [];
end

if nargin < 4 
    filename = false(1);
end

if isfield(opts,'num_frame') num_frame = opts.num_frame; else num_frame = 5; end

idx_frame = round(linspace(1,ntp,num_frame+1));
idx_frame = (idx_frame(1:end-1)+idx_frame(2:end))/2;
t_frame = (idx_frame-1)./(ntp-1);

%% show
im = cell(ntp,1);
% figure('papersize',[7,7],'paperposition',[0,0,7,7]);
% set(gcf,'color','w');

tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
nexttile
for t = 1:ntp
    clf

    plot(x,rho(:,1),'r','LineWidth',2);hold on;
    plot(x,rho1,'b','LineWidth',2);hold on;
    plot(x,rho(:,t),'LineWidth',2);
    pause(0.1);
    drawnow
    frame = getframe(gcf);
    im{t} = frame2im(frame);

end
        

%%
if filename   
    
    % gif
    for idx = 1:ntp
        [A,map] = rgb2ind(im{idx},256);
        if idx == 1
            imwrite(A,map,[filename,'.gif'],'gif','LoopCount',Inf,'DelayTime',0.2);
        else
            imwrite(A,map,[filename,'.gif'],'gif','WriteMode','append','DelayTime',0.2);
        end
    end
    
    
%     % pdf
%     fig = figure('papersize',[5*num_frame,5],'paperposition',[0,0,5*num_frame,5]);
%     figure('papersize',[5*num_frame,5],'paperposition',[0,0,5*num_frame,5]);
%     t = tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
%     nexttile
%     for k = 1:num_frame
%         idx = idx_frame(k);
% 
%         plot(x,rho(:,1),'r','LineWidth',2);hold on;
%         plot(x,rho1,'b','LineWidth',2);hold on;
%         plot(x,rho(:,idx),'LineWidth',2);
% 
%         title(['t=',num2str(t_frame(k))]);
%         exportgraphics(t,[filename,'_',num2str(k),'.pdf'],'BackgroundColor','none')
%     end
%     
%     print(fig,'-dpdf',filename); 


    
end

end