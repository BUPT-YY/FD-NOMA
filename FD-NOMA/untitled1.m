%Simulation for Figure 2.
clc,clear,close all;
npoints = 10;
PUdbm = linspace(10,25,npoints); 
noise = -61.5;%-61.5dBm;
PU = 10.^(PUdbm/10)*10^(-noise/10);

PSI = [100 200 500]; %50, 53, 57dBm;
M = 3; alpha = 4; r1 = 10; r2 = 0;

RU_OMA = zeros(1,length(PUdbm)); Nsim = 100000;
RU_HDNOMA = zeros(1,length(PUdbm));
RU_FDNOMA_500 = zeros(1,length(PUdbm));
RU_FDNOMA_200 = zeros(1,length(PUdbm));
RU_FDNOMA_100 = zeros(1,length(PUdbm));
parfor i=1:length(PU)
    pu = PU(i);
    pr = 0;
    for j=1:Nsim
        x = 0;
        for k=1:M
            hmsqu = untitled0(r1,r2,alpha);
            RU_OMA(i) = RU_OMA(i) + 0.5/M*log2(1+pu*hmsqu);
        end
    end
    for j=1:Nsim
        x = 0;
        for k=1:M
            x = x + untitled0(r1,r2,alpha);
        end
        RU_HDNOMA(i) = RU_HDNOMA(i) + 0.5*log2(1+pu*x);
    end
    for j=1:Nsim
        x = 0;
        for k=1:M
            x = x + untitled0(r1,r2,alpha);
        end
        RU_FDNOMA_100(i) = RU_FDNOMA_100(i) + log2(1+pu*x/(1+PSI(1)*exprnd(1)));
    end
    for j=1:Nsim
        x = 0;
        for k=1:M
            x = x + untitled0(r1,r2,alpha);
        end
        RU_FDNOMA_200(i) = RU_FDNOMA_200(i) + log2(1+pu*x/(1+PSI(2)*exprnd(1)));
    end
    for j=1:Nsim
        x = 0;
        for k=1:M
            x = x + untitled0(r1,r2,alpha);
        end
        RU_FDNOMA_500(i) = RU_FDNOMA_500(i) + log2(1+pu*x/(1+PSI(3)*exprnd(1)));
    end
    RU_OMA(i) = RU_OMA(i)/Nsim;
    RU_HDNOMA(i) = RU_HDNOMA(i)/Nsim;
    RU_FDNOMA_100(i) = RU_FDNOMA_100(i)/Nsim;
    RU_FDNOMA_200(i) = RU_FDNOMA_200(i)/Nsim;
    RU_FDNOMA_500(i) = RU_FDNOMA_500(i)/Nsim;
end

p = plot(PUdbm,RU_OMA,'d-',PUdbm,RU_HDNOMA,'s-',PUdbm,RU_FDNOMA_100,'+-',PUdbm,RU_FDNOMA_200,'o-',PUdbm,RU_FDNOMA_500,'*-','MarkerIndices',linspace(1,npoints,4));
legend('OMA', 'HD-NOMA','FD-NOMA, P_{SI}=100(50dBm)', 'FD-NOMA, P_{SI}=200(53dBm)','FD-NOMA, P_{SI}=500(57dBm)')
xlabel('Transmission Power of Uplink Users(dBm)')
ylabel('Average Uplink Sum Rates')
