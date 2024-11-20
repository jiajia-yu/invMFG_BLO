function gendata(metshape,nt,nx)
    % mesh
    ntp = nt+1; ntm = nt-1;
    dt = 1/nt;

    nxp = nx+1; nxm = nx-1;
    dx = 1/nx;

    x = reshape( linspace(-1/2,1/2,nxp),nxp,1);
    x = (x(1:end-1,:)+x(2:end,:))/2;

    gau = @(mu, sig) exp(-((x-mu)/sig).^2/2)/(sqrt(2*pi)*sig);
    
    switch metshape
        case 'tri'
            g_true = -0.3*cos(2*pi*x) + 0.7;
            lambdaF = 1e-2;
            lambdaG = 5e-1;

            % opts
            opts_gendata = [];
            opts_gendata.maxit = 5e3;
            opts_gendata.tol = 1e-6;
            opts_gendata.stepsize0 = 0.001;
            opts_gendata.stepmodif = 0.8;
            opts_gendata.submaxit = 5;
            opts_gendata.acc = true(1);

            % initial and computation
            N = 2;
            rho0_cell = cell(N,1);
            rho1_cell = cell(N,1);
            rhotilde_cell = cell(N,1); 
            mtilde_cell = cell(N,1);
            rho0_cell{1} = 1.25+0.25*cos(4*pi*(x+0.25));
            rho1_cell{1} = 1.25+0.25*cos(2*pi*x);
            rho0_cell{2} = gau(0,0.1)+0.1;
            rho1_cell{2} = mean(rho0_cell{2})*ones(size(x));
            for n_data = 1:N
                rho0 = rho0_cell{n_data};
                rho1 = rho1_cell{n_data};
%                 tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
%                 nexttile
                plot(x,rho0,'b');hold on; plot(x,rho1,'r'); plot(x,g_true,'k');

                [rho,mx,outs] = solve_mfg_fista(rho0,rho1,g_true,lambdaF,lambdaG,opts_gendata);
                figure;plot(outs.obj_array);
                figure;show_movement(cat(2,rho0,rho),rho1)
                rhotilde_cell{n_data} = rho;
                mtilde_cell{n_data} = mx;
            end
%         case 'tri'
%             g_true = -0.3*cos(2*pi*x) + 0.7;
%             lambdaF = 1e-2;
%             lambdaG = 5e-1;
% 
%             % opts
%             opts_gendata = [];
%             opts_gendata.maxit = 5e3;
%             opts_gendata.tol = 1e-6;
%             opts_gendata.stepsize0 = 0.001;
%             opts_gendata.stepmodif = 0.8;
%             opts_gendata.submaxit = 5;
%             opts_gendata.acc = true(1);
% 
%             % initial and computation
%             N = 1;
%             rho0_cell = cell(N,1);
%             rho1_cell = cell(N,1);
%             rhotilde_cell = cell(N,1); 
%             mtilde_cell = cell(N,1);
%             for n_data = 1:N
%                 rho0 = gau(0,0.1)+0.1;
%                 rho0_cell{n_data} = rho0;
%                 rho1 = mean(rho0)*ones(size(x));
%                 rho1_cell{n_data} = rho1;
% %                 tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
% %                 nexttile
%                 plot(x,rho0,'b');hold on; plot(x,rho1,'r'); plot(x,g_true,'k');
% 
%                 [rho,mx,outs] = solve_mfg_fista(rho0,rho1,g_true,lambdaF,lambdaG,opts_gendata);
%                 figure;plot(outs.obj_array);
%                 figure;show_movement(cat(2,rho0,rho),rho1)
%                 rhotilde_cell{n_data} = rho;
%                 mtilde_cell{n_data} = mx;
%             end
                
        case 'qua'
            g_true = -10*x.^2 + 5;
            lambdaF = 1e-2;
            lambdaG = 5e-1;

            % opts
            opts_gendata = [];
            opts_gendata.maxit = 5e3;
            opts_gendata.tol = 1e-6;
            opts_gendata.stepsize0 = 0.001;
            opts_gendata.stepmodif = 0.8;
            opts_gendata.submaxit = 5;
            opts_gendata.acc = true(1);

            % initial and computation
            N = 2;
            rho0_cell = cell(N,1);
            rho1_cell = cell(N,1);
            rhotilde_cell = cell(N,1); 
            mtilde_cell = cell(N,1);
            
            rho0_cell{1} = gau(0,0.1)+0.1;
            rho1_cell{1} = mean(rho0_cell{1})*ones(size(x));
            rho0_cell{2} = gau(-0.3,0.1) + 0.1;
            rho1_cell{2} = gau( 0.3,0.1) + 0.1;
            for n_data = 1:N
                rho0 = rho0_cell{n_data};
                rho1 = rho1_cell{n_data};
%                 tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
%                 nexttile
                plot(x,rho0,'b');hold on; plot(x,rho1,'r'); plot(x,g_true,'k');

                [rho,mx,outs] = solve_mfg_fista(rho0,rho1,g_true,lambdaF,lambdaG,opts_gendata);
                figure;plot(outs.obj_array);
                figure;
                show_movement(cat(2,rho0,rho),rho1)
                rhotilde_cell{n_data} = rho;
                mtilde_cell{n_data} = mx;
            end
            
        case 'cub'
            g_true = 8*x.*(x-3/8).*(x+3/8) + 1;
            lambdaF = 1e-2;
            lambdaG = 5e-1;

            % opts
            opts_gendata = [];
            opts_gendata.maxit = 5e3;
            opts_gendata.tol = 1e-6;
            opts_gendata.stepsize0 = 0.001;
            opts_gendata.stepmodif = 0.8;
            opts_gendata.submaxit = 5;
            opts_gendata.acc = true(1);

            % initial and computation
            N = 1;
            rho0_cell = cell(N,1);
            rho1_cell = cell(N,1);
            rhotilde_cell = cell(N,1); 
            mtilde_cell = cell(N,1);
            for n_data = 1:N
                rho0 = gau(0,0.1)+0.1;
                rho0_cell{n_data} = rho0;
                rho1 = mean(rho0)*ones(size(x));
                rho1_cell{n_data} = rho1;
%                 tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
%                 nexttile
                plot(x,rho0,'b');hold on; plot(x,rho1,'r'); plot(x,g_true,'k');

                [rho,mx,outs] = solve_mfg_fista(rho0,rho1,g_true,lambdaF,lambdaG,opts_gendata);
                figure;plot(outs.obj_array);
                figure;show_movement(cat(2,rho0,rho),rho1)
                rhotilde_cell{n_data} = rho;
                mtilde_cell{n_data} = mx;
            end
            
        case 'wave'
            g_true = exp(-x-0.5).*sin(4*pi*x)+1.5;
            lambdaF = 1e-2;
            lambdaG = 5e-1;

            % opts
            opts_gendata = [];
            opts_gendata.maxit = 5e3;
            opts_gendata.tol = 1e-6;
            opts_gendata.stepsize0 = 0.001;
            opts_gendata.stepmodif = 0.8;
            opts_gendata.submaxit = 5;
            opts_gendata.acc = true(1);

            % initial and computation
            N = 1;
            rho0_cell = cell(N,1);
            rho1_cell = cell(N,1);
            rhotilde_cell = cell(N,1); 
            mtilde_cell = cell(N,1);
            for n_data = 1:N
                rho0 = gau(0,0.1)+0.1;
                rho0_cell{n_data} = rho0;
                rho1 = mean(rho0)*ones(size(x));
                rho1_cell{n_data} = rho1;
%                 tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
%                 nexttile
                plot(x,rho0,'b');hold on; plot(x,rho1,'r'); plot(x,g_true,'k');

                [rho,mx,outs] = solve_mfg_fista(rho0,rho1,g_true,lambdaF,lambdaG,opts_gendata);
                figure;plot(outs.obj_array);
                figure;show_movement(cat(2,rho0,rho),rho1)
                rhotilde_cell{n_data} = rho;
                mtilde_cell{n_data} = mx;
            end
                
        otherwise
            fprintf('not available');
    end
    save(['../data/',metshape,'_',num2str(nt),num2str(nx),'.mat'],...
         'g_true','lambdaF','lambdaG','rho0_cell','rho1_cell','rhotilde_cell','mtilde_cell');

end