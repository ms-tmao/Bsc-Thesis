% computes the maximal difference between the eigenvalues of a Hermitian
% matrix (such as the obesity matrix Y)
% input: Y Hermitian matrix
% return: real number

function [spread] = spread(Y)
    largest = eigs(Y,1,'largestreal');
    smallest = eigs(Y,1,'smallestreal');
    spread = real(largest-smallest); %the difference should naturally be real; but due to computing inaccuracies it can happen that the imaginary part is not equal to zero
end
                