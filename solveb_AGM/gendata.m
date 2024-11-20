function gendata(obsshape,lambdab,nt,nx,ny)
% mesh
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

gau = @(mu1, mu2, sig1, sig2) exp(-((x-mu1).^2/(sig1^2)+(y-mu2).^2/sig2^2)/2)/(2*pi*sig1*sig2);

switch obsshape
    case 'gaussian'
        b_true = lambdab*gau(0,0,0.1,0.08);
        lambdaF = 1e-1;
        lambdaG = 5e0;

        % opts
        opts_gendata = [];
        opts_gendata.nt = nt;
        opts_gendata.maxit = 5e3;
        opts_gendata.tol = 1e-6;
        opts_gendata.stepsize0 = 0.01;
        opts_gendata.stepmodif = 0.8;
        opts_gendata.submaxit = 5;
        opts_gendata.acc = true(1);

        % initial and computation
        N = 1;
        rho0_cell = cell(N,1);
        rho1_cell = cell(N,1);
        rhotilde_cell = cell(N,1); 
        mtilde_cell = cell(N,2);
        for n_data = 1:N
            rho0 = gau(0, -0.25, 0.08, 0.08) + 0.1;
            rho0_cell{n_data} = rho0;
            rho1 = gau(0,  0.25, 0.08, 0.08) + 0.1;
            rho1_cell{n_data} = rho1;

            [rho,mx,my,outs] = solve_mfg_fista(rho0,rho1,b_true,lambdaF,lambdaG,opts_gendata);
            figure;plot(outs.obj_array);
%                 show_movement(cat(3,rho0,rho),mx,my,b_true/lambdab/10)
            show_movement(cat(3,rho0,rho),mx,my,zeros(size(rho0)))
            rhotilde_cell{n_data} = rho;
            mtilde_cell{n_data,1} = mx;
            mtilde_cell{n_data,2} = my;
        end
        save(['../data/',obsshape,num2str(lambdab,'%.0e'),'_',num2str(nt),num2str(nx),'.mat'],...
            'b_true','lambdaF','lambdaG','rho0_cell','rho1_cell','rhotilde_cell','mtilde_cell');

    case 'twobar'
        % block, terminal and parameters
        b_true = lambdab*(double((x>-0.1) & (x<-0.05) & (y<0)) + double((x>0.05) & (x<0.1) & (y>0)));
        % b_true = double((x>-0.02) & (x<0.02) & (y<-0.2));
        % b_true = double( x.^2+y.^2 < 0.02);
        lambdaF = 1e-1;
        lambdaG = 1e0;

        % opts
        opts_gendata = [];
        opts_gendata.nt = nt;
        opts_gendata.maxit = 5e3;
        opts_gendata.tol = 1e-6;
        opts_gendata.stepsize0 = 0.01;
        opts_gendata.stepmodif = 0.8;
        opts_gendata.submaxit = 5;
        opts_gendata.acc = true(1);

        % initial and computation
        N = 1;
%             initcen = [-0.3,-0.3; -0.3,-0.15; -0.3,0; -0.3,0.15; -0.3,0.3];
        rho0_cell = cell(N,1);
        rho1_cell = cell(N,1);
        rhotilde_cell = cell(N,1); 
        mtilde_cell = cell(N,2);
        for n_data = 1:N
            rho0 = 10*exp(-( (x+0.3).^2 + (y+0.3).^2)*100 ) + 0.1;
            rho0(b_true>0) = 0;
            rho0_cell{n_data} = rho0;
            rho1 = 10*exp(-( (x-0.3).^2 + (y-0.3).^2)*100 ) + 0.1;
            rho1_cell{n_data} = rho1;

            [rho,mx,my,outs] = solve_mfg_fista(rho0,rho1,b_true,lambdaF,lambdaG,opts_gendata);
            figure;plot(outs.obj_array);
            show_movement(cat(3,rho0,rho),mx,my,b_true)
            rhotilde_cell{n_data} = rho;
            mtilde_cell{n_data,1} = mx;
            mtilde_cell{n_data,2} = my;
        end
        save(['../data/',obsshape,num2str(lambdab,'%.0e'),'_',num2str(nt),num2str(nx),'.mat'],...
            'b_true','lambdaF','lambdaG','rho0_cell','rho1_cell','rhotilde_cell','mtilde_cell');

    case 'clover'
        % block, terminal and parameters
        b_true = double( (x-0.19).^2+(y*1.5).^2<0.02 );
        b_true = b_true + double( (x+0.19).^2+(y*1.5).^2<0.02 );
        b_true = b_true + double( (x*1.5).^2+(y-0.19).^2<0.02 );
        b_true = b_true + double( (x*1.5).^2+(y+0.19).^2<0.02 );
        b_true = max(b_true, double( (x).^2+(y).^2<0.005 ));
        b_true = lambdab*b_true;
%             imshow(b_true,[]);colorbar;colormap default
        lambdaF = 1e-1;
        lambdaG = 1e0;

        % opts
        opts_gendata = [];
        opts_gendata.nt = nt;
        opts_gendata.maxit = 5e3;
        opts_gendata.tol = 1e-6;
        opts_gendata.stepsize0 = 0.01;
        opts_gendata.stepmodif = 0.8;
        opts_gendata.submaxit = 5;
        opts_gendata.acc = true(1);

        % initial and computation
        N = 1;
%             initcen = [-0.3,-0.3; -0.3,-0.15; -0.3,0; -0.3,0.15; -0.3,0.3];
        rho0_cell = cell(N,1);
        rho1_cell = cell(N,1);
        rhotilde_cell = cell(N,1); 
        mtilde_cell = cell(N,2);
        for n_data = 1:N
            rho0 = 10*exp(-( (x+0.3).^2 + (y+0.3).^2)*100 ) + 0.1;
            rho0(b_true>0) = 0;
            rho0_cell{n_data} = rho0;
            rho1 = 10*exp(-( (x-0.3).^2 + (y-0.3).^2)*100 ) + 0.1;
            rho1_cell{n_data} = rho1;

            [rho,mx,my,outs] = solve_mfg_fista(rho0,rho1,b_true,lambdaF,lambdaG,opts_gendata);
            figure;plot(outs.obj_array);
%                 show_movement(cat(3,rho0,rho),mx,my,b_true)
            rhotilde_cell{n_data} = rho;
            mtilde_cell{n_data,1} = mx;
            mtilde_cell{n_data,2} = my;
        end
        save(['../data/',obsshape,num2str(lambdab,'%.0e'),'_',num2str(nt),num2str(nx),'.mat'],...
            'b_true','lambdaF','lambdaG','rho0_cell','rho1_cell','rhotilde_cell','mtilde_cell');

    case 'ring'
        radius = sqrt(x.^2+y.^2);
        b_true = double(radius>0.15 & radius<0.25);
        b_true(x>-0.02 & x<0.02) = 0;
        b_true(y>-0.02 & y<0.02) = 0;
        b_true = lambdab*b_true;
%             imshow(b_true,[]);colorbar;colormap default
        lambdaF = 1e-1;
        lambdaG = 5e0;

        % opts
        opts_gendata = [];
        opts_gendata.nt = nt;
        opts_gendata.maxit = 5e3;
        opts_gendata.tol = 1e-6;
        opts_gendata.stepsize0 = 0.01;
        opts_gendata.stepmodif = 0.8;
        opts_gendata.submaxit = 5;
        opts_gendata.acc = true(1);

        % initial and computation
        N = 1;
        rho0_cell = cell(N,1);
        rho1_cell = cell(N,1);
        rhotilde_cell = cell(N,1); 
        mtilde_cell = cell(N,2);
        for n_data = 1:N
            rho0 = gau(-0.3, -0.3, 0.08, 0.08) + 0.1;
            rho0_cell{n_data} = rho0;
            rho1 = gau( 0.3,  0.3, 0.08, 0.08) + 0.1;
            rho1_cell{n_data} = rho1;

            [rho,mx,my,outs] = solve_mfg_fista(rho0,rho1,b_true,lambdaF,lambdaG,opts_gendata);
            figure;plot(outs.obj_array);
            show_movement(cat(3,rho0,rho),mx,my,zeros(size(rho0)))
            rhotilde_cell{n_data} = rho;
            mtilde_cell{n_data,1} = mx;
            mtilde_cell{n_data,2} = my;
        end
        save(['../data/',obsshape,num2str(lambdab,'%.0e'),'_',num2str(nt),num2str(nx),'.mat'],...
            'b_true','lambdaF','lambdaG','rho0_cell','rho1_cell','rhotilde_cell','mtilde_cell');

    otherwise
        fprintf('not available');

end

end