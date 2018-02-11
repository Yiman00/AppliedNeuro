function z = merge(x,y)
% x is a row n-vector with x(1) <= x(2) <= ... <= x(n)
% y is a row m-vector with y(1) <= y(2) <= ... <= y(m) 
% z is a row (m+n)-vector comprised of all the values in x and y, sorted so
% that z(1) <= ... <= z(m+n) 

n = length(x); m = length(y); z = zeros(1,n+m);
ix = 1;     % The index of next x-value to select.
iy = 1;     % The index of next y-value to select. 
for iz = 1:(n+m)
    % Determine the iz-th value for the merged array...
    if ix > n
        % All done with x-values. Select the next y-value.  
        z(iz) = y(iy); iy = iy+1; 
    elseif iy > m
        % All done with y-values. Select the next x-value. 
        z(iz) = x(ix); ix = ix+1; 
    elseif x(ix) <= y(iy) 
        % The next x-value is less than or equal to the next y-value. 
        z(iz) = x(ix); ix = ix+1;
    else
        % The next y-value is less than the next x-value.
        z(iz) = y(iy); iy = iy+1; 
    end
end


        