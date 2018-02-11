function num_moves = hanoicount(n)
% Counts number of moves needed in Towers of Hanoi 
% Example:
%   num_moves = hanoicount(20) for n = 20 pegs

if n==1
    num_moves = 1;
else
    num_moves = 2.*hanoicount(n-1)+1;
end


