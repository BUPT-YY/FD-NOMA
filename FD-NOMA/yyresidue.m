%分子多项式, 分母的根, 分母的各个根的重数
%numerator/denominator = p(x) + sum( fi/gi  )
function [result,p] = yyresidue(numerator, denomRoots, multiplicities)
    L = length(multiplicities);
    [G,g] = polyofRoots(denomRoots, multiplicities);
    f = cell(1,L);
    denominator = conv(G{1,1},g{1,1});
    [p, numerator] = deconv(numerator, denominator);
    if(L==1)
        f{1,1}=numerator; root = denomRoots(1);
        H = horner(numerator, root, multiplicities(1));
        M = multiplicities(1);
        result = zeros(M,3);
        for m=1:M
           result(m,:) = [H(m),root,m];
        end
    else
        for l=1:L
            [~,c2,~] = bezout(g{1,l},G{1,l});
            [~,temp] = deconv(conv(c2,numerator),g{1,l});
            f{1,l} = temp(min(find(temp~=0,1)):end);
        end
        MG = 0;
        for l=1:L
            MG = MG + multiplicities(l);
        end
        result = zeros(MG,3);
        MC = 0;
        for l = 1:L
            M = multiplicities(l);
            root = denomRoots(l);
            H = horner(f{1,l}, root, multiplicities(l));
            for m=1:M
                MC = MC + 1;
                result(MC,:) = [H(m),root,m];
            end
        end
    end
end
%example: 1/(x+1)+3/(x+2)-2/(x+2)^2+4/(x+4)^5+5/(x+10)
%yyresidue([9,256,3118,21180,87444,223676,344296,290128,102560],[-1,-2,-4,-10],[1,2,5,1])
%result:
%   1,-1,1
%   3,-2,1
%  -2,-2,2
%   0,-4,1
%   0,-4,2
%   0,-4,3
%   0,-4,4
%   4,-4,5
%   5,-10,1

%用秦九韶算法(Horner algo.)
%分解 (numerator)/(x-x0)^m = sum( c_j/(x-x0)^j )
function [result] = horner(numerator, denomRoot, denomMultiplicity)
    %补0
    numerator = add(numerator,zeros(1,denomMultiplicity));
    
    if(denomMultiplicity==1)
        result = numerator;
    elseif(denomMultiplicity==2)
        result = [numerator(1),numerator(1)*denomRoot+numerator(2)];
    else
        temp = [1,denomRoot];
        result = [numerator(1)];
        for i=2:denomMultiplicity
            result = conv(result,temp);
            result(end) = result(end)+numerator(i);
        end
    end
end
%例:(13x^2+14x+10)/(x+1)^3 = 13/(x+1) -12/(x+1)^2 + 9/(x+1)^3
%horner([13,14,10],-1,3) -> [13,-12,9]
%例:(-4x^3-45x^2-170x-215)/(x+4)^4 = -4/(x+4)+3/(x+4)^2-2/(x+4)^3+1/(x+4)^4
%horner([-4,-45,-170,-215],-4,4) -> [-4,3,-2,1]


function [G,g] = polyofRoots(denomRoots, multiplicities)
    L = length(multiplicities);
    G = cell(1,L);
    g = cell(1,L);
    for l=1:L
        % (x-x0)^m
        g{1,l} = poly(ones(1, multiplicities(l))*denomRoots(l));
    end
    if (L==1)
        G{1,1} = [1];
    else 
        temp = g{1,1};
        for l=2:L
            temp = conv(temp, g{1,l});
        end
        for l=1:L
            [Gl,~]= deconv(temp, g{1,l});
            G{1,l} = Gl;
        end
    end
    
end
%polyofRoots([-4,-1,-2],[1,3,1])



%一元多项式环中用扩展Euclid算法求Bezout定理的系数
% u(x)f(x)+v(x)g(x) = gcd(f,g)
function [u, v, gcd] = bezout(f, g)
    err = 0.00001;
    index = min(find(f~=0,1)); f = f(index:end);
    index = min(find(g~=0,1)); g = g(index:end);

    um = [1]; u = [0];
    vm = [0]; v = [1];
    if(length(f) < length(g))
        p = f; pm = g;
    else
        p = g; pm = f;
    end
    

    [q,pp] = deconv(pm, p); %一元多项式带余除法 g = qf + r; [q,r] = deconv(g,f)
    gcd = p;
    while(norm(pp)>err)
        up = subtract(um, conv(u,q));
        up = up(min(find(up~=0,1)):end);
        %up = um - conv(u,q);       %长度可能不同, 直接相减有bug
        vp = subtract(vm, conv(v,q));
        vp = vp(min(find(vp~=0,1)):end);
        %vp = vm - conv(v,q);
        um = u; u = up;
        vm = v; v = vp;
        pm = p; p = pp(min(find(pp~=0,1)):end);
    
        gcd = p;
        [q,pp] = deconv(pm, p);
    end
    %转换为首一多项式
    coef = gcd(1);
    gcd = gcd./coef;
    u = u./coef; v = v./coef;
    if(length(f) < length(g))
        temp = u; u = v; v = temp;
    end
    
    
end

function [result] = subtract(a,b)
    m = length(a); n = length(b);
    if(m == n)
        result = a-b;
    elseif(m < n)
        result = -b;
        for i=1:m
            result(n-(m-i)) = result(n-(m-i)) + a(i);
        end
    else
        result = a;
        for i=1:n
            result(m-(n-i)) = result(m-(n-i)) - b(i);
        end
    end
end

function [result] = add(a,b)
    m = length(a); n = length(b);
    if(m == n)
        result = a+b;
    elseif(m < n)
        result = b;
        for i=1:m
            result(n-(m-i)) = result(n-(m-i)) + a(i);
        end
    else
        result = a;
        for i=1:n
            result(m-(n-i)) = result(m-(n-i)) + b(i);
        end
    end
end


%test:
%bezout([1,5,9,7,2],[1,5,6])
