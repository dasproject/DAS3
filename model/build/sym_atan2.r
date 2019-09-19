% sym_atan2.r
% Computes atan2 symbolically
%
% We need this function because Autolev can only do atan2() with numerical arguments
%
% Usage: angle = sym_atan2(y,x)
%
% Method: http://www.mathworks.com/help/toolbox/mupad/stdlib/arg.html
%
atan(#1#/#2#) + pi/2 * sign(#1#) * (1 - sign(#2#))
