% measures of nonnormality
% input: A = a matrix; i in {1,2,3}
% for i=1, measure nu_1 as in the thesis is calculated
% for i=2, measure nu_2 as in the thesis is calculated
% for i=3, departure from normality is calculated
% all defined w.r.t. to Frobenius norm

function [d] = nu(A,i)
    if i==1
        d = norm(A'*A-A*A',"fro");
    end
    if i==2
        E = abs(eig(A));
        S = svd(A);
        d = norm(sort(E,"descend")-S,"inf");
    end
    if i==3
        R = schur(A);
        d = norm(R-diag(diag(R)), "fro");
    end
end
