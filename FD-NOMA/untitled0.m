function [hmsqu]  = untitled0( r1,r2, alpha)
dtoalpha = ((r1*r1-r2*r2)*rand() + r2*r2)^(alpha/2);
hmsqu = exprnd(1)/dtoalpha;