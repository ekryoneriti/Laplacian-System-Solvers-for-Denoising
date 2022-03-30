%forward solve 
function z = forwardSolveBandOrig(L, b, q)
    n = size(L,1);
    z = zeros(n,1);
    for i = 1:n
        z(i) = b(i)/L(i,i);
        for j = max(i-q,1):-1:i
            z(i) = z(i) - L(i,j)*z(j);
        end
    end
end