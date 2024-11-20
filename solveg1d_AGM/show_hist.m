function show_hist(data,opts)

n_im = length(data);
if isfield(opts,'ref') has_ref = true; ref = opts.ref; 
else has_ref = false; end
if isfield(opts,'t_hold') t_hold = opts.t_hold;  
else t_hold = 0.1; end
if isfield(opts,'x_cord') x_cord = opts.x_cord;
else x_cord = linspace(-0.5,0.5,length(data(:,1))); end
if isfield(opts,'it_cord') it_cord = opts.it_cord;
else it_cord = 0:n_im-1; end
if isfield(opts,'filename') savegif = true; filename = opts.filename;
else savegif = false; end


im = cell(n_im,1);
minlim = min(data(:));
maxlim = max(data(:));
if has_ref
    minlim = min(minlim,min(ref(:)));
    maxlim = max(maxlim,max(ref(:)));
end
figure;
set(gcf,'color','w');
for idx = 1:n_im
    clf;
    plot(x_cord,data(:,idx),'r','LineWidth',2);hold on;
    ylim([minlim,maxlim])
    if has_ref
        plot(x_cord,ref,'b','LineWidth',2);
        legend('num','true','location','north');
    end
    title(['number of iteration = ',num2str(it_cord(idx))]);
    pause(t_hold);
    drawnow
    frame = getframe(gcf);
    im{idx} = frame2im(frame);
end

if savegif
    for idx = 1:n_im
        [A,map] = rgb2ind(im{idx},256);
        if idx == 1
            imwrite(A,map,[filename,'.gif'],'gif','LoopCount',Inf,'DelayTime',0);
        else
            imwrite(A,map,[filename,'.gif'],'gif','WriteMode','append','DelayTime',0);
        end
    end
end


end