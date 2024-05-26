% w sweeps of closest_normal.m

function [N, D, R, A_new]=closest_normal_wsweep(A,w)
    n = size(A,1);
    if n < 25
        [N, D, R, A_new]=closest_normal_wsweep_s(A,n,w);
    else
        [N, D, R, A_new]=closest_normal_wsweep_l(A,n,w);
    end
end


function [N, D, R, A_new]=closest_normal_wsweep_s(A,n,w)
    A_new = A;
    R = eye(n,n);
    l=0;
    while l < w
        for j = 1:n
            for k = (j+1):n
                %we save the pivot entries for quick access
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
                A_new = ctranspose(R_new)*A_new*R_new;
                R = R*R_new;
            end
        end
        l=l+1;
    end
    D = diag(diag(A_new));
    N = R * D * ctranspose(R);
end

function [N, D, R, A_new]=closest_normal_wsweep_l(A,n,w)
    A_new = A;
    R = speye(n,n);
    l=0;
    while l < w
        for j = 1:n
            for k = (j+1):n
                %we save the pivot entries for quick access
                a_jj = A_new(j,j);
                a_jk = A_new(j,k);
                a_kj = A_new(k,j);
                a_kk = A_new(k,k);
                % determine shift m and phase psi
                m = 0.5*(a_jj+a_kk);
                pivblockdet = (a_jj-m)*(a_kk-m)-a_kj*a_jk;
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
                A_new = ctranspose(R_new)*A_new*R_new;
                R = R*R_new;
            end
        end
        l=l+1;
    end
    D = diag(diag(A_new));
    N = R * D * ctranspose(R);
end