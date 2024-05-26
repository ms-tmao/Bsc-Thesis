% Our implementation of Ruhe's algorithm
% input: A = a matrix; d = stopping threshold
% return: [N, D, R, A_new, sweep]
% - N is a putative closest normal of A
% - D is diag(A_new) 
% - R is the unitary matrix which transforms A into A_new 
% - A_new is the Delta H-Matrix A is transformed into
% - sweep is the number of sweeps needed for this transformation
%
% We distinguish between small (n<25) and larger matrices. In the second
% case, we use sparse matrices.
function [N, D, R, A_new, sweep]=closest_normal(A, d)
    n = size(A,1);
    if n < 25
        [N, D, R, A_new, sweep]=closest_normal_s(A, d, n);
    else
        [N, D, R, A_new, sweep]=closest_normal_l(A, d, n);
    end
end



function [N, D, R, A_new, sweep]=closest_normal_s(A,d,n)
    R = eye(n,n);
    diff = d + 1;
    sweep = 0;
    while diff > d
        sweep = sweep+1;
        A_new = A;
        for j = 1:n
            for k = (j+1):n
                % save the pivot entries
                a_jj = A_new(j,j);
                a_jk = A_new(j,k);
                a_kj = A_new(k,j);
                a_kk = A_new(k,k);
                % determine shift m and phase psi
                m = 0.5*(a_jj+a_kk);
                pivblockdet = (a_jj-m)*(a_kk-m)-a_kj*a_jk;
                psi = 0.5*angle(-pivblockdet);
                % determine rotation angles phi and alpha
                h_jk = exp(-1i*psi)*a_jk+exp(1i*psi)*conj(a_kj);
                divisor = real(exp(-1i*psi)*(a_jj - a_kk));
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
                A_new = R_new'*A_new*R_new;
                R = R*R_new;
            end
        end
        diff = norm(diag(A_new)) -norm(diag(A));
        A = A_new;
    end
    D = diag(diag(A_new));
    N = R*D*R';
end


function [N, D, R, A_new, sweep]=closest_normal_l(A, d, n)
    R = speye(n,n);
    diff = d + 1;
    sweep = 0;
    while diff > d
        sweep = sweep + 1;
        A_new = A;
        for j = 1:n
            for k = (j+1):n
                % save the pivot entries
                a_jj = A_new(j,j);
                a_jk = A_new(j,k);
                a_kj = A_new(k,j);
                a_kk = A_new(k,k);
                % determine shift m and phase psi
                m = 0.5*(a_jj+a_kk);
                pivblockdet = (a_jj-m)*(a_kk-m)-a_kj*a_jk;
                psi = 0.5*angle(-pivblockdet);
                % determine rotation angles phi and alpha
                h_jk = exp(-1i*psi)*a_jk+exp(1i*psi)*conj(a_kj);
                divisor = real(exp(-1i*psi)*(a_jj - a_kk));
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
                A_new = R_new'*A_new*R_new;
                R = R*R_new;
            end
        end
        diff = norm(diag(A_new))-norm(diag(A));
        A = A_new;
    end
    D = diag(diag(A_new));
    N = R*D*R';
end

