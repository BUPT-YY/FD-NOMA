function [result] = p_eq3(Psi,Pu,r1,r2,alpha,N,M)

    %noise = -61.5;%-61.5dBm;
    %Pu = 10.^(10/10)*10^(-noise/10);Psi = 200;
    
    %N = 20; M = 3; %r1=10;r2=0;alpha = 4;
    b = zeros(1,N); c = zeros(1,N);
    theta = zeros(1,N);
    for i = 1:N
        theta(i) = cos(pi*(2*i-1)/(2*N));
        b(i) = pi/N/(r1+r2)*sqrt(1-theta(i)*theta(i))*((r1+r2)/2 + (r1-r2)/2*theta(i));
        c(i) = ((r1+r2)/2 + (r1-r2)/2*theta(i))^alpha;
    end

    indexs = cell(1, nchoosek(M+N-1,N-1));
    %indexs = zeros(N, nchoosek(M+N-1,N-1));

    stack = [N,M,1];
    while size(stack,1)>0
        top = stack(end,:);
        n = top(1); m = top(2); k = top(3);
        stack(end,:) = [];  %pop
    
        indexs{1,k}(end+1,:)=[N-n+1,m];
        %indexs(N-n+1,k) = m;
        k = k+1;
        if n>1
        for i=1:m
            stack(end+1,:) = [n-1,i,k]; %push
            numsub = nchoosek(i+n-2,n-2);
            for j=1:numsub
                if(m-i~=0)
                    indexs{1,k+j-1}(end+1,:)=[N-n+1,m-i];
                end
                %indexs(N-n+1,k+j-1) = (m-i);
            end
            k = k + numsub;
        end
        end
    end
    
    result = 0;
    
    IS = length(indexs);
    for j=1:IS
        
        %disp(['j= ', num2str(j)]);
        K=indexs{1,j};
        coef = 1; up = M;
        for lk=1:size(K,1)
            coef = coef*nchoosek(up, K(lk,2));
            up = up - K(lk,2);
        end
        numerator = 1; denomRoots = zeros(1,size(K,1));
        multiplicities = zeros(1,size(K,1));
        for lk=1:size(K,1)
            n = K(lk,1); kn = K(lk,2);
            numerator = numerator*((b(n)*c(n))^kn);
            denomRoots(lk) = -c(n);
            multiplicities(lk) = kn;
        end
        [R,~] = yyresidue([numerator], denomRoots,multiplicities);
        temp = 0;
        LR = size(R,1);
        for lr = 1:LR
            A = R(lr,1); cc = -R(lr,2); kmi = R(lr,3);
            if(A~=0)
            tempsub = 0;
            for p=0:(kmi-1)
                tempsub = tempsub + nchoosek(kmi-1,p)*2^(kmi-p)*(phi_eq4(p+kmi+1, cc, Pu, Psi) - phi_eq4(p+kmi, cc, Pu, Psi));
            end
            temp = temp + coef*A*tempsub*exp(-1/Psi)/(Pu^(kmi))/gamma(kmi);
            end
        end
           
        result = result + temp;
        if(mod(j,100)==1)
           disp(['finished ', num2str(100*j/IS), '%']);
        end
    end
end