function [ output_args ] = forwardSolveOrig( input_args )
%BACKSOLVEORIG Summary of this function goes here
%   Detailed explanation goes here
    n = size(L,1);
    z = zeros(n,1);
    for i = 1:n
        z(i) = b(i)/L(i,i);
        for j = i+1:min(i+q,n)
            b(j) = b(j) - L(j,i)*z(i);
        end
    end

end

