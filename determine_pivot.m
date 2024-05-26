% Function to choose the pivot entry minimizing the difference of off(A{s})
% and off(A^{s+1}) at each step s

function [p1, p2, h, divisor, u] = determine_pivot(A, n)
    psi = 0;
    h = 0;
    p1 = 1;
    p2 = 2;
    u = -1;
    for j = 1:n
        for k = j+1:n
            % save the entries
            a_jj = A(j,j);
            a_jk = A(j,k);
            a_kj = A(k,j);
            a_kk = A(k,k);
            % determine shift m and phase psi
            m = 0.5*(a_jj+a_kk);
            d = a_jj-m;
            determinant = -d^2-a_kj*a_jk;
            psi_2 = 0.5*angle(-determinant);
            h_2 = exp(-1i*psi_2)*a_jk+exp(1i*psi_2)*conj(a_kj);
            w_2 = abs(h_2)^2 - 2*imag(exp(-1i*psi_2)*d)^2;
            if w_2 > u
                u = w_2;
                h = h_2;
                psi = psi_2;
                p1 = j;
                p2 = k;
            end
        end
    divisor = real(exp(-1i*psi)*(A(p1,p1) - A(p2,p2)));
    end

     
                