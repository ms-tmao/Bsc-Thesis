% version 2 implementation of Ruhe's algorithm, using a different pivot
% strategy
% input: A = a matrix; z = number of iterations
% return: [N, D, R, A, sweep]
% - N is a putative closest normal of A
% - D is diag(A) 
% - R is the unitary matrix which transforms the input matrix into a Delta H matrix 
% - A is the Delta H-Matrix the input matrix is transformed into
%
% We distinguish between small (n<25) and larger matrices. In the second
% case, we use sparse matrices.
function [N, D, R, A, iter]=closest_normal2(A,d)
    n = size(A,1);
    if n < 25
        [N, D, R, A, iter]=closest_normal2_s(A,n,d);
    else
        [N, D, R, A, iter]=closest_normal2_l(A,n,d);
    end
end



function [N, D, R, A, iter]=closest_normal2_s(A,n,d)
    R = eye(n,n);
    diff = d+1;
    iter = 0;
    while diff > d
        iter = iter+1;
        % choose pivots
        [j, k, h_jk, divisor, u] = determine_pivot(A, n);
        % determine rotation angles phi and alpha
        if h_jk == 0
            phi = 0;
        elseif divisor == 0
            phi = 0.25*pi;
        else
            phi = 0.5*atan(abs(h_jk)/divisor);
        end
        alpha = angle(h_jk);
        % perform rotation
        c = cos(phi);
        s = sin(phi);
        e_value = exp(1i*alpha);
        s1 = e_value*s;
        s2 = conj(e_value)*s;
        R_new = eye(n,n);
        R_new(j,j) = c;
        R_new(j,k) = -s1;
        R_new(k,j) = s2;
        R_new(k,k) = c;
        A = R_new'*A*R_new;
        R = R*R_new;
        diff = u;
    end
    D = diag(diag(A));
    N = R*D*R';
end


function [N, D, R, A, iter]=closest_normal2_l(A,n,d)
    R = speye(n,n);
    diff = d+1;
    iter = 0;
    while diff > d
        iter = iter+1;
        % choose pivots
        [j, k, h_jk, divisor, u] = determine_pivot(A, n);
        % determine rotation angles phi and alpha
        if h_jk == 0
            phi = 0;
        elseif divisor == 0
            phi = 0.25*pi;
        else
            phi = 0.5*atan(abs(h_jk)/divisor);
        end
        alpha = angle(h_jk);
        % perform rotation
        c = cos(phi);
        s = sin(phi);
        e_value = exp(1i*alpha);
        s1 = e_value*s;
        s2 = conj(e_value)*s;
        L = reshape(1:n, n,1);
        L(n+1) = j;
        L(n+2) = k;
        M = reshape(1:n, n,1);
        M(n+1) = k;
        M(n+2) = j;
        V = ones(n,1);
        V(j) = c;
        V(k) = c;
        V(n+1) = -s1;
        V(n+2) = s2;
        R_new = sparse(L,M,V,n,n);
        A = R_new'*A*R_new;
        R = R*R_new;
        diff = u;
    end
    D = diag(diag(A));
    N = R*D*R';
end

