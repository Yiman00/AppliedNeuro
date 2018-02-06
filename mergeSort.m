function y = mergeSort(x)
% x is a vector.
% y is a vector consisting of the values in x sorted from smallest to
% largest. 
%
% Example:
% a = [4 1 6 3 2 9 5 7 6 0];
% b = mergesSort(a); 

n = length(x);
if n==1
    y = x;
else
    m = floor(n/2);
    % Sort the first half...
    y1 = mergeSort(x(1:m)); 
    % Sort the second half...
    y2 = mergeSort(x(m+1:n));
    % Merge...
    y = merge(y1,y2); 
end


    