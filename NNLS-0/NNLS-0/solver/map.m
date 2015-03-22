function [ ax ] = map( x )
global a;

x = reshape(x,625,1);
ax = a*x;

end

