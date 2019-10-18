function W = optimize_informedGL(N,Y, W_multicell, gamma)


[A1,b1,A2,b2,mat_obj] = laplacian_constraint_vech(N);
p = vec(Y*Y')';
norm_x = sum_square_abs(Y);
%% Add the info coming from layers
% sum all the info coming from the multiple views
Layers = sum(cat(3,W_multicell{:}),3);
Layers = Layers + eye(size(Layers));
% Half vectorization
Layers_h = Layers(find(tril(ones(N))));


%% optimization
cvx_begin

% cvx_solver mosek

variable L(N*(N+1)/2,1)

minimize p*mat_obj*L/norm_x + gamma*sum_square_abs(mat_obj*L)

subject to
    A1*L == b1 % zero-row sum and trace constraint
    A2*L <= b2 %non-postive off diag elements constraint
    L(Layers_h == 0) == 0 % constraint on edge set
    %L(non_diag_ind) >= -Layers_max %the weight can be at most the maximum of the layers
cvx_end

%% outputs a vector
L = reshape(mat_obj*L,N,N);
W = -L.*~(eye(size(L))); % get corresponding weight matrix
W(W<10^(-4))=0; % eliminate negligible weights