% PURPOSE
% return KL divergence from P to Q
% 
% INPUT
% P,Q: matrix of column vectors (they must have the same size(,1))
%
% OUTPUT
% kl: matrix of divergences. the first indexes refers to the Q vectors, the
% second regers to the P vectors
%
% REMARK:
% for symmetrized KL use (KL(P,Q)+KL(Q,P))/2

function kl = KLdivergence(P,Q)

dim = size(Q,1);

if (size(P,1) ~= dim)
    error('KLdistance: P and Q dimension mismatch')
end

kl = zeros(size(Q,2),size(P,2));

for i = 1:size(Q,2)
    if mod(i,10) == 0
        fprintf('%d,%d\n',i,size(Q,2))
    end
    index = true(dim,1);
    index(Q(:,i) <= 0) = false;
    this_Q = Q(index,i);
    
    for j = 1:size(P,2)
        this_P = P(index,j);
        kl(i,j) = sum( this_P .* log(this_P./this_Q) );
    end
end

end