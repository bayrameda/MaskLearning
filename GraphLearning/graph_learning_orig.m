function W = graph_learning_orig(N, X, gamma)
% Learning graphs (Laplacian) from structured signals
% Signals X follow Gaussian assumption

max_iter = 30;
alpha = 1;
beta = gamma;

objective = zeros(max_iter,1);
Y_0 = X;
Y = Y_0;
for i = 1:max_iter
    
    % Step 1: given Y, update L
    L = optimize_laplacian_gaussian(N,Y,alpha,beta);
        
    % Step 2: Given L, update Y
    R = chol(eye(N) + alpha*L);
    Y = R \ (R' \ (Y_0));
    
    % plot the objective
    % objective(i) = norm(Y-Y_0,'fro')^2 + alpha*trace(Y'*L*Y) + beta*(norm(L,'fro')^2);
    objective(i) = norm(Y-Y_0,'fro')^2 + alpha*vec(Y*Y')'*vec(L) + beta*(norm(L,'fro')^2);
%     figure(3)
%     plot(i,objective(i), '.r');
%     hold on, drawnow
    
    % stopping criteria
    if i>=2 && abs(objective(i)-objective(i-1))<10^(-4)
        display(['Learning stopped at iteration-', num2str(i)]);
        break
    end
    
end
% Convert
W = -L.*~(eye(size(L))); % get corresponding weight matrix
W(W<10^(-4))=0; % eliminate negligible weights