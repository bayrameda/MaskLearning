clc, clear all, close all;
addpath(genpath('../'));
load('Meteo.mat')
num_edge= @(W) sum(W(:)>0)/2;
set(0,'defaulttextInterpreter','latex');
cvx_setup;
cvx_solver Mosek
%% Measurements
meas_indx = 1; % choose measurement type from records_name
Observation = cell2mat(cellfun(@(x) x(meas_indx,:),Meteo.records_81_10(:),'UniformOutput',false));
st_indx = find(~all(Observation'==0));
year_avg = Observation(st_indx,13);
X = Observation(st_indx,1:12);
N = numel(st_indx); % number of stations
Meteo.st_indx = st_indx;
Wmulti = generateLayers(Meteo);
%% ML
[M_multi, W_mle, W_m, W_e] = optimize_MLextended(N, X, Wmulti, 10^7);
A_gps = M_multi{1};
A_alt = M_multi{2};
display(['Measurement type:', Meteo.records_name(meas_indx)]);
gps_percentage = num_edge(A_gps)/ num_edge(W_mle)
alt_percentage = num_edge(A_alt)/ num_edge(W_mle)
%% Reorganize wrt year average
[~, sort_ind]=sort(year_avg);
W_ord_gps = Wmulti{1}(sort_ind, sort_ind);
W_ord_alt = Wmulti{2}(sort_ind, sort_ind);
W_ord_gpsMask = M_multi{1}(sort_ind, sort_ind);
W_ord_altMask = M_multi{2}(sort_ind, sort_ind);
% Sparsity patterns
figure;
set(gcf, 'Position', [1725 541 405 161]);
subplot(121);
% GPS in blue
s = imagesc(2*(W_ord_gps>0));
s.AlphaData = (W_ord_gps>0);
hold on;
% Alt in red
s = imagesc(3*(W_ord_alt>0));
s.AlphaData = 0.5*(W_ord_alt>0);
s.AlphaData((W_ord_alt>0) > (W_ord_gps>0)) = 1;
colormap hsv; set(gca,'YTick',[], 'XTick',[]);
title('Layers');
subplot(122);
% GPS in blue
s = imagesc(2*(W_ord_gpsMask>0));
s.AlphaData = (W_ord_gpsMask>0);
hold on;
% Alt in red
s = imagesc(3*(W_ord_altMask>0));
s.AlphaData = 0.5*(W_ord_altMask>0);
s.AlphaData((W_ord_altMask>0) > (W_ord_gpsMask>0)) = 1;
colormap hsv; set(gca,'YTick',[], 'XTick',[]);
title('Inferred Masks');
hold on;
h = zeros(2, 1);
h(1) = bar(NaN,'b');
h(2) = bar(NaN,'r');
lgn = legend(h, 'GPS','Altitude',...
   'Orientation','horizontal','Location','southoutside'); 
lgn.Interpreter = 'latex';