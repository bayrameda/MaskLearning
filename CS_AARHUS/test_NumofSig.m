clc, clear all, close all;
addpath(genpath('../'));
cvx_setup
set(0,'defaulttextInterpreter','latex');
load CS_lunch_fb_work.mat
num_edge= @(W) sum(W(:)>0)/2;
precision = @(W,Wgt) sum(W(:) & Wgt(:))/sum(W(:)>0);
recall = @(W,Wgt) sum(W(:) & Wgt(:))/sum(Wgt(:)>0);
F1 = @(p,r) 2*p*r/(p+r);
%%
W_lunch = CSmult.lunch;
W_fb = CSmult.fb;
W_work = CSmult.work;
Coord = CSmult.coord;
N = length(W_lunch);
W_multicell = {W_fb,W_work};
%% change number of signals
X_all = CSmult.signal;
num_sig = size(X_all,2);
n_sig = round(linspace(2,num_sig, 10));
n_rep = [15*ones(9,1);1]; %repeat each exp for 15 times
gamma = [repmat(10^7,1,5), repmat(0.6,1,5)];
vol = sum(W_lunch(:));
f_gli =zeros(size(n_sig));
f_glorig = zeros(size(n_sig));
f_ml =zeros(size(n_sig));
%% At each loop learning on different number of signals
for ii= 1:numel(n_sig)
    glorig = [];
    gli = [];
    mle = [];
    for jj = 1:n_rep(ii)
        idx_sig = randsample(num_sig,n_sig(ii));
        X = X_all(:,idx_sig);
        %% GL-SigRep
        betalfa = 0.015;
        W_gl = graph_learning_orig(N, X, betalfa);
        p = precision(W_gl,W_lunch);
        r = recall(W_gl, W_lunch);
        f = F1(p,r);
        glorig = [glorig ; p, r, f];
        %% GL-informed
        betalfa = 0.012;
        W_gli = optimize_informedGL(N, X, W_multicell, betalfa);
        p = precision(W_gli,W_lunch);
        r = recall(W_gli, W_lunch);
        f = F1(p,r);
        gli = [gli ; p, r, f];
        %% ML
        [ext_M_multi, W_mle, W_m, W_e] = optimize_MLextended(N, X, W_multicell, gamma(ii),vol);
        p = precision(W_mle,W_lunch);
        r = recall(W_mle, W_lunch);
        f = F1(p,r);
        mle = [mle ; p, r, f];
    end
    f_glorig(ii) = mean(glorig(:,3));
    f_gli(ii) = mean(gli(:,3));
    f_ml(ii) = mean(mle(:,3));
end
%% Draw
figure;
plot(n_sig, f_ml, 'b.-');
hold on;
plot(n_sig, f_gli, 'r.-');
plot(n_sig, f_glorig, 'k.-');
ylabel('F-score');
xlabel('Number of signals');
title('F-score vs Number of Signals');
legend('ML','GL-informed', 'GL-SigRep');