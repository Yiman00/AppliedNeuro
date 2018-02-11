% Fibonacci Term using Recursion

function [result] = fibonacci_term(n)

if n==0||n==1
    result = n;
else
    result = fibonacci_term(n-2)+fibonacci_term(n-1);
end
end

