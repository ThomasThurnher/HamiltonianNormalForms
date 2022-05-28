function [Hk] = linear_composition_coefficients(A,H,k)
% Calculates the composition coefficients of power series as given in 
% https://doi.org/10.1016/j.jsv.2020.115640
% Appendix C

l = size(H{1},2);
z_k = nchoosek(k+l-1,l-1);
N = size(H{1},1);
K = flip(sortrows(nsumk(l,k,'nonnegative')).',2);

string = 'revlex';
Hk = zeros(N,z_k,k);

%find nonzero minima of all multi-indices of order k
[kj,Ikj] = min(K +(k+1)*(K ==0)); 
%%
% Loop over all orders of multi-indices that are potentially nonzero and do
% not require the knowledge of the order $k$ coefficients of the SSM parametrisation.
for ord = 1:1
    %%
    % Create all multi-indices of order $k-1$
    if l >1
        g       = flip(sortrows(nsumk(l,k-1,'nonnegative')).',2); %lexicographic!!
    else 
        g = k-1;
    end
    
    [gmink,I_k,I_g] = multi_subtraction(K,g,'Parametrised');
    
    %%
    % Now all the positions of the order $m$ multi-indices are calculated. Assume
    % the vector $\texttt{gmink}$ has $u_\alpha$ columns.
    I_gmk     = multi_index_2_ordering(gmink,string,[]);
    %%
    % We read out all the values $u_j$ positions of minimal entries
    pos     = Ikj(I_k);                 %find row for every vector in mv
    sz      = size(I_k,2);
    poslin  = sub2ind([l,sz],pos,1:sz); %create linear index array
    % read out minimal entries
    gmk_j     = gmink(poslin);            %has to be a row vector
    %%
    % One thing that is redundant is all the calculations where $u_j$ is zero. The
    % corresponding values could be ignored, by taking the indices corresponding to
    % $u_j$= 0 out of calculations .
    %
    % Consequently, the coefficients $u_jS_{i,\mathbf{u}} H_{i,s-1,\mathbf{k}_f-\mathbf{u}}
    % $  are calculated. In the array $\texttt{Hk}$ the first row corresponds to the
    % phase space dimension, the third one to the scalar exponents.
    %
    % $\texttt{H}_{k-u}(:,\texttt{I}_{g},:)$ is a matrix in $\mathbb{C}^{2n \times
    % \texttt{sz} \times k-1}$. Contribution of coefficients where $s=1$ explicitly
    % depend on the SSM coefficients of the order that is calculated and will be added
    % in the end after the highest order coefficients of the SSM parametrisation have
    % been computed,but fortunately are not needed for the calculations of themselves.
    % Finally $\texttt{S\{ord\}}(:,\texttt{gmk\_I}) \in \mathbb{C}^{2n \times \texttt{sz}
    % }$.
    
    %read out contributions
    coeff   = A(:,I_gmk).* gmk_j .* H{k-ord}(:,I_g,:);

    
    %%
    % $\texttt{coeff} \in \mathbb{C}^{2n \times \texttt{sz} \times k-1}$. What is
    % left now is to add up all the parts corresponding to a specific multi-index
    % of order  $k$ and multiply with the $s$ the slice corresponds to. We want to
    % add all the columns corresponding to where $\texttt{I}_{k}$ has the same values
    % (the columns contain coefficients that contribute to the same order $k$ multi-index).
    for f = unique(I_k)
        Hk(:,f,2:(k-ord+1)) = Hk(:,f,2:(k-ord+1)) + sum(coeff(:,I_k==f,:),2);
    end
end
%multiplication with s and kj
s  = 1:k;
%match dimensions for pointiwse multiplication
Hk = permute(Hk,[2,3,1]) .* ((1./kj).' * s);
Hk = permute(Hk,[3,1,2]);

end


