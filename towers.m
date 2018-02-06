function [] = towers(n, frompeg, topeg, auxpeg)
% Towers of Hanoi 
% Example:
%   towers(5, 'A','C','B') where n = 5 pegs 

if n==1
    fprintf('\t move disk 1 from peg %c to peg %c \n', frompeg, topeg);
else
    towers(n-1,frompeg,auxpeg,topeg);
    fprintf('\t move disk %d from peg %c to peg %c \n', n, frompeg, topeg);
    towers(n-1,auxpeg,topeg,frompeg);
end




