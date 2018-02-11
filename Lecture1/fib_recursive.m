% Fibonacci Sequence (Recursion)

% Provide base case
n(1) = 1;
n(2) = 1;

k = 3;
while k <= 10
    n(k) = n(k-1)+n(k-2);
    k = k+1;
end


