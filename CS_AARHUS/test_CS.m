clc, clear all, close all;
addpath(genpath('../'));
cvx_setup
set(0,'defaulttextInterpreter','latex');
load CS_lunch_fb_work.mat
num_edge= @(W) sum(W(:)>0)/2;
precision = @(W,Wgt) sum(W(:) & Wgt(:))/sum(W(:)>0);
recall = @(W,Wgt) sum(W(:) & Wgt(:))/sum(Wgt(:)>0);
jaccard = @(W1,W2) sum(W1(:)&W2(:))/sum(W1(:)|W2(:));
F1 = @(p,r) 2*p*r/(p+r);
%%
W_lunch = CSmult.lunch;
W_fb = CSmult.fb;
W_work = CSmult.work;
Coord = CSmult.coord;
X = CSmult.signal;
N = length(W_lunch);
W_multicell = {W_fb,W_work};
%% Draw
figure;
subplot(131);
plot(graph(W_fb),'b','XData', Coord(:,1),'YData',Coord(:,2),'NodeLabel',{});
xlim([0,1]), ylim([0,1]);
title(['FB, $\left|\mathcal{E}\right|$ =', num2str(sum(W_fb(:)>0)/2)]);
subplot(132);
plot(graph(W_work),'r','XData', Coord(:,1),'YData',Coord(:,2),'NodeLabel',{});
xlim([0,1]), ylim([0,1]);
title(['Work, $\left|\mathcal{E}\right|$ =', num2str(sum(W_work(:)>0)/2)]);
subplot(133);
plot(graph(W_lunch),'XData', Coord(:,1),'YData',Coord(:,2),'NodeLabel',{});
xlim([0,1]), ylim([0,1]);
title(['Lunch, $\left|\mathcal{E}\right|$ =', num2str(sum(W_lunch(:)>0)/2)]);
%% Stats
union = num_edge(W_fb + W_work)
global_not_in_union = num_edge(W_lunch > (W_fb | W_work))
common = num_edge((W_fb & W_work))
common_not_in_global = num_edge((W_fb & W_work) > W_lunch)
%% GL-SigRep
betalfa = 0.015;
W_gl = graph_learning_orig(N, X, betalfa);
p_G = precision(W_gl,W_lunch);
r_G = recall(W_gl, W_lunch);
f_G_gl = F1(p_G,r_G);
glorig_graph = [p_G, r_G, f_G_gl];
figure;
plot(graph(W_gl),'XData', Coord(:,1),'YData',Coord(:,2),'NodeLabel',{});
xlim([0,1]), ylim([0,1]);
title(['GL-SigRep, $\left|\mathcal{E}\right|$ =', num2str(num_edge(W_gl)), ' precision-recall-F1 = ', num2str(glorig_graph)]);
display('/////////////GL-SigRep is done!/////////////////');
%% GL-informed
betalfa = 0.012;
W_gli = optimize_informedGL(N, X, W_multicell, betalfa);
p_gli = precision(W_gli,W_lunch);
r_gli = recall(W_gli, W_lunch);
f_gli = F1(p_gli,r_gli);
gli = [p_gli, r_gli, f_gli];
figure;
plot(graph(W_gli),'XData', Coord(:,1),'YData',Coord(:,2),'NodeLabel',{});
xlim([0,1]), ylim([0,1]);
title(['GL-informed, $\left|\mathcal{E}\right|$ =', num2str(num_edge(W_gli)), ' precision-recall-F1 = ', num2str(gli)]);
%% ML
gamma = 0.6;
vol = sum(W_lunch(:));
[ext_M_multi, W_mle, W_m, W_e] = optimize_MLextended(N, X, W_multicell, gamma,vol);
p_mle = precision(W_mle,W_lunch);
r_mle = recall(W_mle, W_lunch);
f_mle = F1(p_mle,r_mle);
mle = [p_mle,r_mle, f_mle];
% Draw ML result
figure;
set(gcf, 'Position', [228 506 1235 237]);
subplot(141);
plot(graph(W_mle & ext_M_multi{1}),'-ob','XData', Coord(:,1),'YData',Coord(:,2),'NodeLabel',{});
xlim([0,1]), ylim([0,1]);
e_fb = num_edge(W_mle & ext_M_multi{1});
title(['FB, $\left|\mathcal{E}\right|$ =', num2str(e_fb)]);
subplot(142);
plot(graph(W_mle & ext_M_multi{2}),'-or','XData', Coord(:,1),'YData',Coord(:,2),'NodeLabel',{});
xlim([0,1]), ylim([0,1]);
e_work = num2str(num_edge(W_mle & ext_M_multi{2}));
title(['Work, $\left|\mathcal{E}\right|$ =', num2str(e_work)]);
subplot(143)
plot(graph((W_e>0) > (W_m>0)),'-og','XData', Coord(:,1),'YData',Coord(:,2),'NodeLabel',{});
xlim([0,1]), ylim([0,1]);
e_add = num_edge((W_e>0) > (W_m>0));
title(['Additive Correction, $\left|\mathcal{E}\right|$ =', num2str(e_add)]);
subplot(144);
plot(graph(W_mle),'XData', Coord(:,1),'YData',Coord(:,2),'NodeLabel',{});
xlim([0,1]), ylim([0,1]);
title(['MLext, $\gamma$ = ', num2str(gamma), ', $\left|\mathcal{E}\right|$ =', num2str(num_edge(W_mle)),...
    ', precision-recall-F1 = ', num2str(mle)]);
%% ML global vs Ground-truth global
figure;
subplot(121);
plot(graph(W_lunch&W_fb),'-ob', 'XData', Coord(:,1),'YData',Coord(:,2),'NodeLabel',{});
hold on;
plot(graph(W_lunch&W_work),'-or', 'XData', Coord(:,1),'YData',Coord(:,2),'NodeLabel',{});
plot(graph(W_lunch>(W_work+W_fb)),'-og', 'XData', Coord(:,1),'YData',Coord(:,2),'NodeLabel',{});
xlim([0,1]), ylim([0,1]);
title(['Ground-truth Lunch, $\mathcal{E}_{Lunch}\cap \mathcal{E}_{FB}$=', num2str(num_edge(W_lunch&W_fb)),...
    ', $\mathcal{E}_{Lunch}\cap \mathcal{E}_{Work}$=', num2str(num_edge(W_lunch&W_work)), ', additional=', num2str(global_not_in_union)]);
subplot(122);
plot(graph(W_mle & ext_M_multi{1}),'-ob','XData', Coord(:,1),'YData',Coord(:,2),'NodeLabel',{});
hold on;
plot(graph(W_mle & ext_M_multi{2}),'-or','XData', Coord(:,1),'YData',Coord(:,2),'NodeLabel',{});
plot(graph((W_e>0) > (W_m>0)),'-og','XData', Coord(:,1),'YData',Coord(:,2),'NodeLabel',{});
xlim([0,1]), ylim([0,1]);
title(['MLext, mask-FB=', num2str(e_fb), ', mask-Work=', num2str(e_work), ', additional=', num2str(e_add)]);
%% Sparsity patterns
figure;
subplot(121);
% fb in blue
s = imagesc(2*(W_lunch&W_fb));
s.AlphaData = (W_lunch&W_fb);
hold on;
% work in red
s = imagesc(3*(W_lunch&W_work));
s.AlphaData = 0.5*(W_lunch&W_work);
s.AlphaData((W_lunch&W_work) > (W_lunch&W_fb)) = 1;
% additional in green
s = imagesc((W_lunch>(W_work+W_fb)));
s.AlphaData = (W_lunch>(W_work+W_fb));
colormap hsv;
title('Ground-truth lunch');
subplot(122);
% fb in blue
s = imagesc(2*(W_mle & ext_M_multi{1}));
s.AlphaData = (W_mle & ext_M_multi{1});
hold on;
% work in red
s = imagesc(3*(W_mle & ext_M_multi{2}));
s.AlphaData = 0.5*(W_mle & ext_M_multi{2}); % mixed color on the common
s.AlphaData ((W_mle & ext_M_multi{2}) > (W_mle & ext_M_multi{1})) = 1;
% additional in green
s = imagesc((W_e>0) > (W_m>0));
s.AlphaData = (W_e>0) > (W_m>0);
colormap hsv;
title('ML-ext lunch');
%% Groundtruth vs ML result
figure;
subplot(221);
plot(graph(W_lunch&W_fb),'-ob', 'XData', Coord(:,1),'YData',Coord(:,2),'NodeLabel',{});
hold on;
plot(graph(W_lunch&W_work),'-or', 'XData', Coord(:,1),'YData',Coord(:,2),'NodeLabel',{});
plot(graph(W_lunch>(W_work+W_fb)),'-og', 'XData', Coord(:,1),'YData',Coord(:,2),'NodeLabel',{},'NodeColor','k');
xlim([0,1]), ylim([0,1]); set(gca,'YTick',[], 'XTick',[]);
title(['Ground truth lunch graph, $\left|\mathcal{E}\right|$ = ', num2str(num_edge(W_lunch))]);
subplot(222);
s = imagesc(2*(W_lunch&W_fb));
s.AlphaData = (W_lunch&W_fb);
hold on;
% work in red
s = imagesc(3*(W_lunch&W_work));
s.AlphaData = 0.5*(W_lunch&W_work);
s.AlphaData((W_lunch&W_work) > (W_lunch&W_fb)) = 1;
% additional in green
s = imagesc((W_lunch>(W_work+W_fb)));
s.AlphaData = (W_lunch>(W_work+W_fb));
colormap hsv; set(gca,'YTick',[], 'XTick',[]);
title('Ground-truth lunch');
hold on;
h = zeros(4, 1);
h(1) = bar(NaN,'b');
h(2) = bar(NaN,'r');
h(3) = bar(NaN,'FaceColor',[0.4940 0.1840 0.5560]);
h(4) = bar(NaN,'g'); 
lgn = legend(h, ['Facebook, $\left|\mathcal{E}\right|$ = ', num2str(num_edge(W_lunch&W_fb))],...
    ['Work, $\left|\mathcal{E}\right|$ = ', num2str(num_edge(W_lunch&W_work))],...
    'Common',['Additional, $\left|\mathcal{E}\right|$ = ', num2str(num_edge(W_lunch>(W_work+W_fb)))],...
    'Orientation','horizontal','Location','southoutside'); lgn.Interpreter = 'latex';
%%%%%%%% ML
subplot(223);
plot(graph(W_mle & ext_M_multi{1}),'-ob','XData', Coord(:,1),'YData',Coord(:,2),'NodeLabel',{});
hold on;
plot(graph(W_mle & ext_M_multi{2}),'-or','XData', Coord(:,1),'YData',Coord(:,2),'NodeLabel',{});
plot(graph((W_e>0) > (W_m>0)),'-og','XData', Coord(:,1),'YData',Coord(:,2),'NodeLabel',{},'NodeColor','k');
xlim([0,1]), ylim([0,1]); set(gca,'YTick',[], 'XTick',[]);
title(['ML, $\left|\mathcal{E}\right|$ =', num2str(num_edge(W_mle))]);
subplot(224);
s = imagesc(2*(W_mle & ext_M_multi{1}));
s.AlphaData = (W_mle & ext_M_multi{1});
hold on;
% work in red
s = imagesc(3*(W_mle & ext_M_multi{2}));
s.AlphaData = 0.5*(W_mle & ext_M_multi{2}); % mixed color on the common
s.AlphaData ((W_mle & ext_M_multi{2}) > (W_mle & ext_M_multi{1})) = 1;
% additional in green
s = imagesc((W_e>0) > (W_m>0));
s.AlphaData = (W_e>0) > (W_m>0);
colormap hsv; set(gca,'YTick',[], 'XTick',[]);
title('ML sparsity pattern');
hold on;
h = zeros(4, 1);
h(1) = bar(NaN,'b');
h(2) = bar(NaN,'r');
h(3) = bar(NaN,'FaceColor',[0.4940 0.1840 0.5560]);
h(4) = bar(NaN,'g'); 
lgn = legend(h, ['Facebook, $\left|\mathcal{E}\right|$ = ', num2str(num_edge(ext_M_multi{1}&W_mle))],...
    ['Work, $\left|\mathcal{E}\right|$ = ', num2str(num_edge(ext_M_multi{2}&W_mle))],...
    'Common',['Additional, $\left|\mathcal{E}\right|$ = ', num2str(num_edge((W_e>0) > (W_m>0)))],...
    'Orientation', 'horizontal','Location','southoutside'); lgn.Interpreter = 'latex';





