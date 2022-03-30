%forward solve 
function x = backwardSolveBand(U, z, q)
n = size(U,1);
x = zeros(n,1);
for i = n:-1:1
    x(i) = z(i)/U(i,i);
    for j = max(i-q,1):i-1
        z(j) = z(j) - U(j,i)*x(i);
    end
end
end