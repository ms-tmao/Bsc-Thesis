% w iterations of closest_normal2.m

function [N, D, R, A]=closest_normal2_witer(A,w)
    n = size(A,1);
    if n < 25
        [N, D, R, A]=closest_normal2_witer_s(A,n,w);
    else
        [N, D, R, A]=closest_normal2_witer_l(A,n,w);
    end
end



function [N, D, R, A]=closest_normal2_witer_s(A,n,w)
    R = eye(n,n);
    for z = 1:w
        % choose pivots
        [j, k, h_jk, divisor] = determine_pivot(A, n);
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
    end
    D = diag(diag(A));
    N = R*D*R';
end


function [N, D, R, A]=closest_normal2_witer_l(A,n,w)
    R = speye(n,n);
    for z = 1:w
        % choose pivots
        [j, k, h_jk, divisor] = determine_pivot(A, n);
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
    end
    D = diag(diag(A));
    N = R*D*R';
end

