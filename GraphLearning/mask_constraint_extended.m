function [dup_W, dup_L  mat_obj1, mat_obj2, row_sum, off_diag, trace] = mask_constraint_extended(N,T)

%% matrix for objective (vech -> vec)
dup_W = DuplicationM_d(N, 'lo_d'); %indices to duplicate the elements
% (N^2, N*(N-1)/2) mapping of indices to the lu indices
% fills the main diagonal with zeros
dup_L = DuplicationM_d(N, 'lo'); %together with diagonal
%% matrix to convert full weight vector to full degree vector
I = vec(repmat(linspace(1,N*N,N),N,1));
J = 1:N*N;
D_dup = sparse(I,J,ones(N*N,1));

%% matrix to convert half vectorized conc. weight matrices to full vectorized conc. weight matrices
Dup_cell = cell(T,1);
Dup_cell(:) = {dup_W};
mat_obj1 = blkdiag(Dup_cell{:});

%% matrix to convert full vectorized conc. weight matrices to full vectorized conc. degree matrices
D_dup_cell = cell(T,1);
D_dup_cell(:) = {D_dup};
mat_obj2 = blkdiag(D_dup_cell{:});

%% Row sum matrix to be multiplied with (N X N+1 /2, 1)
X = ones(N);
[r,c] = size(X);
i     = 1:numel(X);
j     = repmat(1:c,r,1);
B     = sparse(i',j(:),X(:))';
row_sum = B*dup_L;

%% non-positive off-diagonal entries to be multiplied with (N X N+1 /2, 1)
for i = 1:N
    tmp{i} = ones(1,N+1-i);
    tmp{i}(1) = 0;
end
off_diag = spdiags(horzcat(tmp{:})',0,N*(N+1)/2,N*(N+1)/2); %off-diags

%% Trace to be multiplied with (N X N+1 /2, 1)
trace = sparse(ones(1,N*(N+1)/2)-horzcat(tmp{:}));
