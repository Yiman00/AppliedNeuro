% Compute the cross variance between x and y
% 
% Input
% x = signal #1
% y = signal #2 
%
% Output
% rxy = cross variance between x and y
% lag = lag axis

function [rxy, lag] = cross_variance(x,y)

N = length(x); % size of data
rxyP = zeros(1, N-1); % cross-variance at + shifts
rxpN = zeros(1, N-1); % cross-varance at - shifts

lagP = zeros(1, N-1); % lag axis for + shifts
lagN = zeros(1, N-1); % lag axis for - shifts

for h = 1:N-1 % + shifts
    temp = sum( circshift(x, [1,h]).*y); 
    rxyP(h) = temp;
    lagP(h) = h;
end

for h = 1:N-1 % - shifts
    temp = sum( circshift(x, [1,-h]).*y); 
    rxyN(N-1-(h-1)) = temp;
    lagN(N-1-(h-1)) = -h; 
end

temp = sum(x.*y); % zero shift
rxy0 = temp;

rxy = [rxyN rxy0 rxpP]/N; % final output
lag = [lagN 0 lagP];
end







    




