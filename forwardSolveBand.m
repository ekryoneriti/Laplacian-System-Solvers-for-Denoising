%forward solve 
function z = forwardSolveBand(L, b, q)
n = size(L,1);
z = zeros(n,1);
for i = 1:n
    z(i) = b(i)/L(i,i);
    for j = i+1:min(i+q,n)
        z(j) = z(j) - L(j,i)*z(i);
    end
end
end