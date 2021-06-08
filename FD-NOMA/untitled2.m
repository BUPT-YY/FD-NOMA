%Simulation for Figure 3.
clc,clear,close all;

npoints = 4;
PUdbm = linspace(10,25,npoints);
noise = -61.5;%-61.5dBm;
PU = 10.^(PUdbm/10)*10^(-noise/10);

PSI = 200;


PrU1 = zeros(1,length(PU));%M=2,r1=10,r2=0
PrU2 = zeros(1,length(PU));%M=3,r1=10,r2=0
PrU3 = zeros(1,length(PU));%M=2,r1=15,r2=5
PrU4 = zeros(1,length(PU));%M=3,r1=15,r2=5

M = 2;r1=10;r2=0;alpha=4;
Nsim = 1000000;
parfor k=1:length(PU)
    pu=PU(k);
for i = 1:Nsim
    x = 0;
    for m=1:M
        x = x + untitled0(r1,r2,alpha);
    end
    y = exprnd(1);
    if pu*x+1<= PSI*PSI*y*y
        PrU1(k) = PrU1(k) + 1;
    end
end
    PrU1(k) = PrU1(k)/Nsim;
end

M = 3;r1=10;r2=0;alpha=4;
parfor k=1:length(PU)
    pu=PU(k);
for i = 1:Nsim
    x = 0;
    for m=1:M
        x = x + untitled0(r1,r2,alpha);
    end
    y = exprnd(1);
    if pu*x+1<= PSI*PSI*y*y
        PrU2(k) = PrU2(k) + 1;
    end
end
    PrU2(k) = PrU2(k)/Nsim;
end

M = 2;r1=15;r2=5;alpha=4;
parfor k=1:length(PU)
    pu=PU(k);
for i = 1:Nsim
    x = 0;
    for m=1:M
        x = x + untitled0(r1,r2,alpha);
    end
    y = exprnd(1);
    if pu*x+1<= PSI*PSI*y*y
        PrU3(k) = PrU3(k) + 1;
    end
end
    PrU3(k) = PrU3(k)/Nsim;
end

M = 3;
parfor k=1:length(PU)
    pu=PU(k);
for i = 1:Nsim
    x = 0;
    for m=1:M
        x = x + untitled0(r1,r2,alpha);
    end
    y = exprnd(1);
    if pu*x+1<= PSI*PSI*y*y
        PrU4(k) = PrU4(k) + 1;
    end
end
    PrU4(k) = PrU4(k)/Nsim;
end

PrU1_theory = zeros(1,length(PU));%M=2,r1=10,r2=0
PrU2_theory = zeros(1,length(PU));%M=3,r1=10,r2=0
PrU3_theory = zeros(1,length(PU));%M=2,r1=15,r2=5
PrU4_theory = zeros(1,length(PU));%M=3,r1=15,r2=5
N=20;
parfor k=1:length(PU)
    pu=PU(k);
    PrU1_theory(k) = p_eq3(PSI,pu,10,0,4,N,2);
    PrU2_theory(k) = p_eq3(PSI,pu,10,0,4,N,3);
    PrU3_theory(k) = p_eq3(PSI,pu,15,5,4,N,2);
    PrU4_theory(k) = p_eq3(PSI,pu,15,5,4,N,3);
end


semilogy(PUdbm, PrU1,'d-',PUdbm,PrU2,'*-',PUdbm, PrU3,'o-',PUdbm,PrU4,'s-',PUdbm,PrU1_theory,'-.',PUdbm,PrU2_theory,'-.',PUdbm,PrU3_theory,'-.',PUdbm,PrU4_theory,'-.','MarkerIndices',linspace(1,npoints,4))
legend('M = 2, r_{1} = 10, r_2 = 0','M = 3, r_{1} = 10, r_2 = 0','M = 2, r_{1} = 15, r_2 = 5','M = 3, r_{1} = 15, r_2 = 5','Location','SouthWest')
xlabel('Transmission Power of Uplink Users(dBm)')
ylabel('Probability for HD outperforming FD')
