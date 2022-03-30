function x = columnOrientedForSolve(L, b)
    bandwidth = size(L,1);
    x = zeros(size(b));
    for i = 1:bandwidth
        x(i) = b(i)/L(i,i);
        for j = i+1:bandwidth
            b(j) = b(j) - L(j,i)*x(i);
        end
    end
end