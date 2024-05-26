% We use this function to create Jordan matrices 
% input: array ns (contains sizes of the Jordan blocs), array diagvals
% (contains eigenvalues of Jordan blocs; ns and diagvals must have same size)
% returns: Jordan matrix with size(ns) blocs of eigenvalue as in diagvals

% N = size of the whole matrix, for N > 25 we again use sparse matrices 

function [J, U] = create_Jordan_matrix(ns, diagvals)
    N = sum(ns);
    D = diagvals(1)*ones(1,ns(1));
    k = 1;
    while k < size(ns,2)
        k = k+1;
        D = [D diagvals(k)*ones(1,ns(k))];
    end
    U = ones(1,N-1);
    k=0;
    i=1;
    while i < size(ns,2)
        k = k+ns(1,i);
        U(1,k) = 0;
        i=i+1;
    end
    if N < 26
        U2 = [diag(U); zeros(1,N-1)];
        U2 = [zeros(N,1) U2];
        J = diag(D)+U2;
    else
        J = sparse(N,N);
        L = [reshape(1:N,N,1); reshape(1:N-1,N-1,1)];
        M = [reshape(1:N,N,1); reshape(2:N,N-1,1)];
        V = [D U];
        J = sparse(L,M,V, N, N);
    end
end


