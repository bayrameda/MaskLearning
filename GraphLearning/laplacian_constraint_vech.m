function [A1,b1,A2,b2,mat_obj] = laplacian_constraint_vech(N)
% A1 : matrix to give row-sum and diag-sum when multiplied with Laplacian
% matrix
% b1 : equality values A1*L =b1
% A2 : matrix to give the off diagonal elements of Laplacian matrix
% b2 : inequality values A2*L <= b2
% mat_obj : matrix to map the lower diagonal indices (half vector) to the vector indices
% mat_obj*L converts vech(Lmat) to vec(Lmat)
%
% constraints        
% mat_cons1*L == zeros(N,1)
% mat_cons2*L <= 0
% vec_cons3*L == N


%% matrix for objective (vech -> vec)
mat_obj = DuplicationM(N); %indices to duplicate the elements
% (N^2, N*(N+1)/2) mapping of indices to the lu indices

%% matrix for constraint 1 (zero row-sum)
% B = sparse(N,N^2);
% for i = 1:N
%     B(i,(i-1)*N+1:i*N) = ones(1,N);
% end
X = ones(N);
[r,c] = size(X);
i     = 1:numel(X);
j     = repmat(1:c,r,1);
B     = sparse(i',j(:),X(:))';
mat_cons1 = B*mat_obj; 

%% matrix for constraint 2 (non-positive off-diagonal entries)
for i = 1:N
    tmp{i} = ones(1,N+1-i);
    tmp{i}(1) = 0;
end
mat_cons2 = spdiags(horzcat(tmp{:})',0,N*(N+1)/2,N*(N+1)/2); %off-diags

%% vector for constraint 3 (trace constraint)
vec_cons3 = sparse(ones(1,N*(N+1)/2)-horzcat(tmp{:}));

%% create constraint matrices
% equality constraint A2*vech(L)==b2
% row-sum should be zero and diag sum should equal to N
A1 = [mat_cons1;vec_cons3]; %indices for row-sum and diag-sum
b1 = [sparse(N,1);N]; % zeros and trace


% inequality constraint A1*vech(L)<=b1
% off diags should be non-positive
A2 = mat_cons2; %indices 
b2 = sparse(N*(N+1)/2,1); % zeros