function M = DuplicationM_d(n, option)
% function M = DuplicationM(n)
% M = DuplicationM(n, 'lo') % (default) OR
% M = DuplicationM(n, 'up') % 
% M = DuplicationM(n, 'lo_d') % 
% Return duplication matrix order n
%
%
% Ouput are sparse
%
% DuplicationM(size(A),'lo')*A(itril(size(A))) == A(:) %true for lower half
% DuplicationM(size(A),'up')*A(itriu(size(A))) == A(:) %true for upper half
% DuplicationM(size(A),'lo_d')*A(itril(size(A),-1)) == A(:) %true for
% lowerhalf excluding the main diagonal, and fill it with zeros.
%
% Author: Bruno Luong <brunoluong@yahoo.com>
% Date: 21/March/2009
%
% Ref : Magnus, Jan R.; Neudecker, Heinz (1980), "The elimination matrix:
% some lemmas and applications", Society for Industrial and Applied Mathematics.
% Journal on Algebraic and Discrete Methods 1 (4): 422–449,  
% doi:10.1137/0601049, ISSN 0196-5212.

if nargin<2
    option = 'lo'; % default
end

if isscalar(n)
    n = [n n];
end

switch option
    case 'lo' 
        [I J] = itril(n);
    case 'up' 
        [I J] = itriu(n);
    case 'lo_d' 
        [I J] = itril(n, -1);
    otherwise
        error('option must be ''lo'' or ''up''');
end

% Find the sub/sup diagonal part that can flip to other side
loctri = find(I~=J & J<=n(1) & I<=n(2));
% Indices of the flipped part
Itransposed = sub2ind(n, J(loctri), I(loctri));

% Convert to linear indice
I =  sub2ind(n, I, J);

% Result
M = sparse([I; Itransposed], ...
           [(1:length(I))'; loctri], 1, prod(n), length(I));

end
