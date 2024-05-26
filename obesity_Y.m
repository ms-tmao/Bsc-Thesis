% computes the obesity matrix Y to a Delta H-matrix as in definition 2
% input: A (must  have Delta H-matrix form)
% return: obesity matrix Y of A

function [Y] = obesity_Y(A) %A is in H matrix form
    n = size(A,1);
    Y = zeros(n,n);
    for j=1:n
        for k=1:n
            if j~=k && A(k,k)-A(j,j) ~= 0
                Y(j,k) = A(j,k)/(A(k,k)-A(j,j));
            end
        end
    end
                