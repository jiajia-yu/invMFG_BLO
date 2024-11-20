% load('results\wave_0e+00_N=1_T=5000_K=10.mat')
% load('results\cub_1e-05_N=1_T=5000_K=10.mat')
% load('results\qua_0e+00_N=1_T=5000_K=10.mat')
% load('results\qua_0e+00_N=2_T=5000_K=10.mat')
% load('results\qua_1e-04_N=1_T=5000_K=10.mat')
load('results\qua_1e-04_N=2_T=5000_K=10.mat')



opts_show.ref = g_true;
opts_show.t_hold = 0.001;
opts_show.x_cord = x;
opts_show.it_cord = 0:10:size(g_hist,2);
opts_show.filename = ['videos\',filename];
show_hist(g_hist(:,1:10:end),opts_show);
