clear
clc

filename = 'gaussian1e+00_1664_new_N=1_T=6000_K=5.mat';
% filename = 'gaussian5e-01_1664_new_N=1_T=6000_K=5.mat';
% filename = 'gaussian1e-01_1664_new_N=1_T=6000_K=5.mat';
% filename = 'gaussian5e-02_1664_new_N=1_T=6000_K=5.mat';

load(['results\',filename])
opts_gendata.maxit = 2000;
%%
for count_t = 1:length(b_hist)
    rho_tmp = cell(N,1);
    m_tmp = cell(N,2);
    for n_data = 1:N
        [rho_tmp{n_data},m_tmp{n_data,1},m_tmp{n_data,2},outs] = ...
            solve_mfg_fista(rho0_cell{n_data},rho1_cell{n_data},b_hist{count_t},lambdaF,lambdaG,opts_gendata);
    end
    valD_hist(count_t) = comp_valD(dd,rho_tmp,m_tmp,rhotilde_cell,mtilde_cell);
end
save(['results/',filename,'_new_N=',num2str(N),'_T=',num2str(T),'_K=',num2str(K)],...
      'valD_hist','-append');

