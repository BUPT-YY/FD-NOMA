

function [ result ] = paraCylinder(nu, z)
    if (nu > 0)
        result = 2^(nu/2)*exp(-z*z/4)/sqrt(pi)*(cos(pi*nu/2)*gamma((nu+1)/2)*hypergeom([-nu/2],[1/2],z*z/2)+sqrt(2)*z*sin(pi*nu/2)*gamma(nu/2+1)*hypergeom([1/2-nu/2],[3/2],z*z/2));
    elseif(nu == 0)
        result = exp(-z*z/4);
    else
        k = -nu - 1;
        result = 2^((k-1)/2)/gamma(k+1)*exp(-z*z/4)*(gamma((k+1)/2)*hypergeom([(k+1)/2],[1/2],z*z/2)-sqrt(2)*z*gamma(k/2+1)*hypergeom([k/2+1],[3/2],z*z/2));
    end
end