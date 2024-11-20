function gendata(nt,nx,ny)

%% mesh
ntp = nt+1; ntm = nt-1;
dt = 1/nt;

nxp = nx+1; nxm = nx-1;
dx = 1/nx;

nyp = ny+1; nym = ny-1;
dy = 1/ny;

x = reshape( linspace(-1/2,1/2,nxp),nxp,1);
x = repmat(x, 1,nyp);
x = (x(1:end-1,:)+x(2:end,:))/2;
x = (x(:,1:end-1)+x(:,2:end))/2;
y = reshape( linspace(-1/2,1/2,nyp),1,nyp);
y = repmat(y, nxp,1);
y = (y(1:end-1,:)+y(2:end,:))/2;
y = (y(:,1:end-1)+y(:,2:end))/2;

%% forward solver
metshape = '2dtri';
g_true = (cos(2*pi*(x-0.25)).*sin(2*pi*(y))+1.5)/2;
lambdaF = 1e-1;
lambdaG = 1e0;

% opts
opts_gendata = [];
opts_gendata.maxit = 5e3;
opts_gendata.tol = 1e-6;
opts_gendata.stepsize0 = 0.001;
opts_gendata.stepmodif = 0.8;
opts_gendata.submaxit = 5;
opts_gendata.acc = true(1);

% initial and computation
N = 5;
initcen = [0.3,-0.3; 0, -0.3; -0.3,-0.3; -0.3,0; -0.3,0.3];
termcen = [-0.3,0.3; 0,  0.3;  0.3, 0.3;  0.3,0; 0.3,-0.3];
rho0_cell = cell(N,1);
rho1_cell = cell(N,1);
rhotilde_cell = cell(N,1); 
mtilde_cell = cell(N,2);
for n_data = 1:N
    rho0 = 10*exp(-( (x-initcen(n_data,1)).^2 + (y-initcen(n_data,2)).^2)*100 ) + 0.1;
    rho1 = 10*exp(-( (x-0.3).^2 + (y-0.3).^2)*100 ) + 0.1;
%     rho1 = 10*exp(-( (x-termcen(n_data,1)).^2 + (y-termcen(n_data,2)).^2)*100 ) + 0.1;
    tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
    nexttile
    imshow(max(rho0,rho1),[]);colormap default
    rho0_cell{n_data} = rho0;
    rho1_cell{n_data} = rho1;

    [rho,mx,my,outs] = solve_mfg_fista(rho0,rho1,{g_true+4,g_true+2,g_true+1},lambdaF,lambdaG,opts_gendata);
    figure;plot(outs.obj_array);
    show_movement(cat(3,rho0,rho),mx,my,[],['illu/',metshape,'_',num2str(n_data)])
    rhotilde_cell{n_data} = rho;
    mtilde_cell{n_data,1} = mx;
    mtilde_cell{n_data,2} = my;
end
save(['../data/',metshape,'_',num2str(nt),num2str(nx),'.mat'],...
    'g_true','lambdaF','lambdaG','rho0_cell','rho1_cell','rhotilde_cell','mtilde_cell');

end