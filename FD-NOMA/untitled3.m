%Simulation for Figure 4.
%clc,clear,close all


noise = -61.5;
pu = 10.^(20/10)*10^(-noise/10);%PU = 20dBm, Noise = -61.5dBm

npoints = 5;%93
PBSdbm = linspace(20,40,npoints);
PBS = 10.^(PBSdbm./10)*10^(-noise/10);

d0 = 100;Nsim = 100000;

result0 = zeros(1,npoints);
M = 2; r1 = 10; r2 = 0; alpha = 4; theta = pi/8;
parfor k=1:npoints
    pbs = PBS(k);
    for i=1:Nsim
        z = 0; x = exprnd(1/d0^alpha);
        for j=1:M
            z = z + exprnd(1/d0^alpha);
        end
        if x*pbs+1 <= pu*pu*z*z
            result0(k) = result0(k) + 1;
        end
    end
    result0(k) = result0(k)/Nsim;
end

result0ana = zeros(1,npoints);
parfor k=1:npoints
   pbs = PBS(k);
   result0ana(k) = 1-lowerIncompleteGamma(M,d0^alpha/pu)/gamma(M);
   for p=0:M-1
       result0ana(k) = result0ana(k) - d0^(alpha*M)*exp(-d0^alpha/pu)/gamma(M)/(pu^M)*nchoosek(M-1,p)*phi_eq4(p+1, d0^alpha, pbs, pu/(d0^alpha));
   end
end



result1 = zeros(1,npoints);
M = 3; r1 = 10; r2 = 0; alpha = 4;
parfor k=1:npoints
    pbs = PBS(k);
    for i=1:Nsim
        z = 0; x = exprnd(1/d0^alpha);
        for j=1:M
            z = z + exprnd(1/d0^alpha);
        end
        if x*pbs+1 <= pu*pu*z*z
            result1(k) = result1(k) + 1;
        end
    end
    result1(k) = result1(k)/Nsim;
end
result1ana = zeros(1,npoints);
parfor k=1:npoints
   pbs = PBS(k);
   result1ana(k) = 1-lowerIncompleteGamma(M,d0^alpha/pu)/gamma(M);
   for p=0:M-1
       result1ana(k) = result1ana(k) - d0^(alpha*M)*exp(-d0^alpha/pu)/gamma(M)/(pu^M)*nchoosek(M-1,p)*phi_eq4(p+1, d0^alpha, pbs, pu/(d0^alpha));
   end
end


%
result2 = zeros(1,npoints);
M = 2; r1 = 10; r2 = 0; alpha = 4;
parfor k=1:npoints
    pbs = PBS(k);
    for i=1:Nsim
        z = 0; x = exprnd(1/d0^alpha);
        for j=1:M
            r = r2 + (r1-r2)*rand();
            phi = (theta/2)*rand();
            d = sqrt((d0+r*cos(phi))^2+(r*sin(phi))^2);
            z = z + exprnd(1/d^alpha);
        end
        if x*pbs+1 <= pu*pu*z*z
            result2(k) = result2(k) + 1;
        end
    end
    result2(k) = result2(k)/Nsim;
end

%
result3 = zeros(1,npoints);
M = 3; r1 = 10; r2 = 0; alpha = 4;
parfor k=1:npoints
    pbs = PBS(k);
    for i=1:Nsim
        z = 0; x = exprnd(1/d0^alpha);
        for j=1:M
            r = r2 + (r1-r2)*rand();
            phi = (theta/2)*rand();
            d = sqrt((d0+r*cos(phi))^2+(r*sin(phi))^2);
            z = z + exprnd(1/d^alpha);
        end
        if x*pbs+1 <= pu*pu*z*z
            result3(k) = result3(k) + 1;
        end
    end
    result3(k) = result3(k)/Nsim;
end

%
result4 = zeros(1,npoints);
M = 2; r1 = 15; r2 = 5; alpha = 4;
parfor k=1:npoints
    pbs = PBS(k);
    for i=1:Nsim
        z = 0; x = exprnd(1/d0^alpha);
        for j=1:M
            r = r2 + (r1-r2)*rand();
            phi = (theta/2)*rand();
            d = sqrt((d0+r*cos(phi))^2+(r*sin(phi))^2);
            z = z + exprnd(1/d^alpha);
        end
        if x*pbs+1 <= pu*pu*z*z
            result4(k) = result4(k) + 1;
        end
    end
    result4(k) = result4(k)/Nsim;
end

%
result5 = zeros(1,npoints);
M = 3; r1 = 15; r2 = 5; alpha = 4;
parfor k=1:npoints
    pbs = PBS(k);
    for i=1:Nsim
        z = 0; x = exprnd(1/d0^alpha);
        for j=1:M
            r = r2 + (r1-r2)*rand();
            phi = (theta/2)*rand();
            d = sqrt((d0+r*cos(phi))^2+(r*sin(phi))^2);
            z = z + exprnd(1/d^alpha);
        end
        if x*pbs+1 <= pu*pu*z*z
            result5(k) = result5(k) + 1;
        end
    end
    result5(k) = result5(k)/Nsim;
end

semilogy(PBSdbm, result0,'x-',PBSdbm, result0ana,'x-.',PBSdbm, result1,'o-',PBSdbm, result1ana,'x-.',PBSdbm, result2,'+-',PBSdbm, result3,'s-',PBSdbm, result4,'d-',PBSdbm, result5,'*-','MarkerIndices',linspace(1,npoints,5));
legend('Upper bound, M=2','Upper bound, M=2, ana','Upper bound, M=3','Upper bound, M=3, ana','M = 2, r_{1} = 10, r_2 = 0','M = 3, r_{1} = 10, r_2 = 0','M = 2, r_{1} = 15, r_2 = 5','M = 3, r_{1} = 15, r_2 = 5','Location','SouthWest')
xlabel('Transmission Power of the BS(dBm)')
ylabel('Probability for HD outperforming FD')


function [result] = lowerIncompleteGamma(s,x)
    result = gamma(s) - igamma(s,x);
end