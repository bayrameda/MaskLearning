function [M_multicell, W_global, W_m, W_e] = optimize_MLextended(N, Y, W_multicell, gamma, varargin)
% Mask solver for MLextended
T = numel(W_multicell);
getHalfVec = @(Mat) Mat(find(tril(ones(N),-1)));
W_veccell = cellfun(getHalfVec, W_multicell, 'UniformOutput',false);
Wmulti = vec(cell2mat(W_veccell));
vol = N;
if nargin == 5
    vol = varargin{1};
end
%% Laplacian constraints
[dup_W, dup_L  mat_obj1, mat_obj2, row_sum, off_diag, trace] = mask_constraint_extended(N,T);
p = vec(Y*Y')';
p_multi = repmat(p,1,T);

d_m = N*(N-1)*T/2; % dimension of vectorized masks
d_e = N*(N+1)/2; % dimension of vectorized laplacian_additional
B_m = [eye(d_m), zeros(d_m,d_e)]; %(d_m X d_m+d_e) draw mask variable
B_e = [zeros(d_e,d_m), eye(d_e)]; %(d_e X d_m+d_e) draw laplacian_additional variable
l2w = diag(off_diag) > 0 ;
% S = [M;L_e]

%% optimization

if (gamma <= 10^6) %relaxed
    cvx_begin

    variable S(d_m + d_e,1)

    minimize  p_multi*(mat_obj2*mat_obj1*(Wmulti.*B_m * S) - mat_obj1*(Wmulti.*B_m * S)) + p*(dup_L*B_e * S) + gamma*sum_square_abs(dup_L*B_e * S)% sum tr(XLtX) + tr(XL_eX) + frob(L_e)

    subject to
        2*Wmulti'*B_m * S + trace*B_e * S == vol %  volume constraint
        B_m * S >= 0 %non-negative masks
        sum(reshape(B_m * S,N*(N-1)/2,T),2) == 1 %edge_loc %bound constraint for layer mask sum
        row_sum*B_e * S == 0 %zero row sum of laplacian additional
        %off_diag*B_e * S <= 0 %non-positive off-diagonal for laplacian additional
        sum(reshape(Wmulti.*B_m * S,N*(N-1)/2,T),2) - B_e(l2w,:)*S >= 0 % non-negative element for global weight matrix


    cvx_end
    
    % check solved
    assert(abs(2*Wmulti'*B_m * S + trace*B_e * S - vol) < 10^(-4), 'Not solved!');

    S(abs(S) < 10^(-4)) = 0;
    M = B_m * S; %half-vector of masks
    L_e = reshape(dup_L*B_e * S,N,N);
    
else %Constrained problem
    
    cvx_begin
    
    variable M(d_m,1)

    minimize  p_multi*(mat_obj2*mat_obj1*(Wmulti.*M) - mat_obj1*(Wmulti.*M))% sum tr(XLtX)

    subject to
        2*Wmulti'*M == vol %  volume constraint
        M >= 0 %non-negative masks
        sum(reshape(M,N*(N-1)/2,T),2) == 1 %edge_loc %bound constraint for layer mask sum

    cvx_end
    %check solved
    assert(abs(2*Wmulti'*M -vol) < 10^(-4), 'Not solved!');
    M(M < 10^(-4)) = 0;
    L_e = zeros(N);

end

M(Wmulti < 10^(-5)) = 0; % zero out the mask associated non-existing edges
%% convert global weight half vector to matrix form
W_m = reshape(dup_W*sum(reshape(Wmulti.*M,N*(N-1)/2,T),2),N,N);

W_e = -L_e.*~(eye(N)); % get corresponding weight matrix for additional

W_global = W_m + W_e;

W_global(W_global < 10^(-4)) = 0;

% outputs a vector
%% convert mask half vectors 
M_multicell = num2cell(reshape(mat_obj1*M, N*N,T),1);
shapeN = @(Mat) reshape(Mat,N,N);
M_multicell = cellfun(shapeN, M_multicell,'UniformOutput',false);


