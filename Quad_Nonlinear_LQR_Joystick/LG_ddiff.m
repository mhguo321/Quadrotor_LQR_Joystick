function [ x2 ] = LG_ddiff( v,Fs )
%LG_DDIFF Summary of this function goes here
%   Detailed explanation goes here
% [ x2 ] = LG_ddiff( v,Fs )
% Î¢·ÖÆ÷
%

h = 1/Fs;
L = length(v);
r = 20.0;
x1 = v;
x2 = zeros(L,1);
for i=2:L
    x1(i) = x1(i-1) + h * x2(i-1);
    x2(i) = x2(i-1) + h *( -(r^2)*(x1(i-1) - v(i-1)) - 2*r*x2(i-1) );
end
end

