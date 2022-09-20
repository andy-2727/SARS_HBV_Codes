clear
clc
tic
% Latin Hypercube Sampling (LHS)/ PRCC implementation for HPV_TB Co-infection model
global   Lambda mu2 beta_C2 beta_H2 beta_CH2 varthet_12 varthet_22 alpha_12 alpha_22 eta_C2 eta_H2 eta_CH2 zet_12 zet_22 xi_C2 xi_H2 xi_CH2 
Lambda = 732500/(76.31*365); % this parameter is fixed
N = 1000; %number of points or samples or runs
p = 16; % number of parameters to be sampled.
%      1       2    3      4      5     6     7     8      9    10    11    12    13     14     15     16               
lb = [0.00001 0.1   0.1   0.1   0.1   0.1    0.01  0.01  0.01  0.01  0.01   0.5   0.5    0.03   0.03   0.03]; %lower bounds of the 5 parameters
ub = [0.00005 2.0   2.0   2.0   0.8   0.8    1.0   1.0   0.05  0.05  0.05   1.5   1.5    0.5    0.5    0.5]; %upper bound of the 5 parameters
X = lhsdesign(N,p,'criterion','correlation'); % generate uniformly distributed samples in a normalized design
D = bsxfun(@plus,lb,bsxfun(@times,X,(ub-lb))); %maps the result of X into the interval determine by ub and lb
%calculate output response. In this case, it is the reproduction number
% However, we can use the complete set of parameter values and simulate the
% system of ODE and make use of any of the state variables as an output
% response.
%Lets start, using the reproduction number. For the model of interest (see
%below), this is
Ro = D(:,2).*(D(:,5).*(D(:,9) + D(:,14) + D(:,1))+D(:,7))./((D(:,7)+D(:,1)).*(D(:,9) + D(:,14) + D(:,1))) %SARS-CoV-2
%Ro = D(:,3).*(D(:,6).*(D(:,10) + D(:,15) + D(:,1))+D(:,8))./((D(:,8)+D(:,1)).*(D(:,10) + D(:,15) + D(:,1))) % HBV
%Ro = D(:,4)./(D(:,11) + D(:,16) + D(:,1)) % co-infection
% %We now proceed with performing the sensitivity analysis using the PRCC
% %technique.
% % We first rank the LHS matrix, D, and the output matrix (in this case a
% % vector, Ro).
 [Ds,Di] = sort(D);
 [Dx,Dr] = sort(Di);
% %Dr is now the rank-transformed matrix of D
% %Next we rank Ro
[Ros, Roi] = sort(Ro);
[Rox, Ror] = sort(Roi);
% %Ror is now the rank-transformed vector of Ro
 Dr;
 Ror;
% % We keep track of the respective parameters of the model from the matrix D
mu = Dr(:,1);
beta_C = Dr(:,2);
beta_H = Dr(:,3);
beta_CH = Dr(:,4);
varthet_1 = Dr(:,5);
varthet_2 = Dr(:,6);
alpha_1 = Dr(:,7);
alpha_2 = Dr(:,8);
eta_C = Dr(:,9);
eta_H = Dr(:,10);
eta_CH = Dr(:,11);
zet_1 = Dr(:,12);
zet_2 = Dr(:,13);
xi_C = Dr(:,14);
xi_H = Dr(:,15);
xi_CH = Dr(:,16);
% 
% % We carry out two linear regressions per parameter and outcome and obtain
% % the residuals from the regressions
% % For beta
Xmu = [ones(size(mu)) beta_C beta_H  varthet_1 eta_C alpha_1 beta_CH alpha_2 eta_H xi_H xi_C zet_2 zet_1 varthet_2 eta_CH];
[muReg,bintmu,rmu] = regress(mu,Xmu);
[RoRegmu,bintRomu,rRomu] = regress(Ror,Xmu);
% For beta_C
Xbeta_C = [ones(size(beta_C)) beta_H  mu varthet_1 eta_C alpha_1 beta_CH alpha_2 eta_H xi_H xi_C zet_2 zet_1 varthet_2 eta_CH];
[beta_CReg,bintbeta_C,rbeta_C] = regress(beta_C,Xbeta_C);
[RoRegbeta_C,bintRobeta_C,rRobeta_C] = regress(Ror,Xbeta_C);
% For beta_H
Xbeta_H = [ones(size(beta_H)) beta_C  mu varthet_1 eta_C alpha_1 beta_CH alpha_2 eta_H xi_H xi_C zet_2 zet_1 varthet_2 eta_CH];
[beta_HReg,bintbeta_H,rbeta_H] = regress(beta_H,Xbeta_H);
[RoRegbeta_H,bintRobeta_H,rRobeta_H] = regress(Ror,Xbeta_H);
% For beta_CH
Xbeta_CH = [ones(size(beta_CH)) beta_C beta_H  mu varthet_1 eta_C alpha_1 alpha_2 eta_H xi_H xi_C zet_2 zet_1 varthet_2 eta_CH];
[beta_CHReg,bintbeta_CH,rbeta_CH] = regress(beta_CH,Xbeta_CH);
[RoRegbeta_CH,bintRobeta_CH,rRobeta_CH] = regress(Ror,Xbeta_CH);
% For varthet_1
Xvarthet_1 = [ones(size(varthet_1)) beta_C beta_H  mu eta_C alpha_1 beta_CH alpha_2 eta_H xi_H xi_C zet_2 zet_1 varthet_2 eta_CH];
[varthet_1Reg,bintvarthet_1,rvarthet_1] = regress(varthet_1,Xvarthet_1);
[RoRegvarthet_1,bintRovarthet_1,rRovarthet_1] = regress(Ror,Xvarthet_1);
% For varthet_2
Xvarthet_2 = [ones(size(varthet_2)) beta_C beta_H  mu varthet_1 eta_C alpha_1 beta_CH alpha_2 eta_H xi_H xi_C zet_2 zet_1 eta_CH];
[varthet_2Reg,bintvarthet_2,rvarthet_2] = regress(varthet_2,Xvarthet_2);
[RoRegvarthet_2,bintRovarthet_2,rRovarthet_2] = regress(Ror,Xvarthet_2);
% For alpha_1
Xalpha_1 = [ones(size(alpha_1)) beta_C beta_H  mu varthet_1 eta_C beta_CH alpha_2 eta_H xi_H xi_C zet_2 zet_1 varthet_2 eta_CH];
[alpha_1Reg,bintalpha_1,ralpha_1] = regress(alpha_1,Xalpha_1);
[RoRegalpha_1,bintRoalpha_1,rRoalpha_1] = regress(Ror,Xalpha_1);
% For alpha_2
Xalpha_2 = [ones(size(alpha_2)) beta_C beta_H  mu varthet_1 eta_C alpha_1 beta_CH eta_H xi_H xi_C zet_2 zet_1 varthet_2 eta_CH];
[alpha_2Reg,bintalpha_2,ralpha_2] = regress(alpha_2,Xalpha_2);
[RoRegalpha_2,bintRoalpha_2,rRoalpha_2] = regress(Ror,Xalpha_2);
% For eta_C
Xeta_C = [ones(size(eta_C)) beta_C beta_H  mu varthet_1 alpha_1 beta_CH alpha_2 eta_H xi_H xi_C zet_2 zet_1 varthet_2 eta_CH];
[eta_CReg,binteta_C,reta_C] = regress(eta_C,Xeta_C);
[RoRegeta_C,bintRoeta_C,rRoeta_C] = regress(Ror,Xeta_C);
% For eta_H
Xeta_H = [ones(size(eta_H)) beta_C beta_H  mu varthet_1 eta_C alpha_1 beta_CH alpha_2 xi_H xi_C zet_2 zet_1 varthet_2 eta_CH];
[eta_HReg,binteta_H,reta_H] = regress(eta_H,Xeta_H);
[RoRegeta_H,bintRoeta_H,rRoeta_H] = regress(Ror,Xeta_H);
% For eta_CH
Xeta_CH = [ones(size(eta_CH)) beta_C beta_H  mu varthet_1 eta_C alpha_1 beta_CH alpha_2 eta_H xi_H xi_C zet_2 zet_1 varthet_2];
[eta_CHReg,binteta_CH,reta_CH] = regress(eta_CH,Xeta_CH);
[RoRegeta_CH,bintRoeta_CH,rRoeta_CH] = regress(Ror,Xeta_CH);
% For zet_1
Xzet_1 = [ones(size(zet_1)) beta_C beta_H  mu varthet_1 eta_C alpha_1 beta_CH alpha_2 eta_H xi_H xi_C zet_2 varthet_2 eta_CH];
[zet_1Reg,bintzet_1,rzet_1] = regress(zet_1,Xzet_1);
[RoRegzet_1,bintRozet_1,rRozet_1] = regress(Ror,Xzet_1);
% For zet_2
Xzet_2 = [ones(size(zet_2)) beta_C beta_H  mu varthet_1 eta_C alpha_1 beta_CH alpha_2 eta_H xi_H xi_C zet_1 varthet_2 eta_CH];
[zet_2Reg,bintzet_2,rzet_2] = regress(zet_2,Xzet_2);
[RoRegzet_2,bintRozet_2,rRozet_2] = regress(Ror,Xzet_2);
% For xi_C
Xxi_C = [ones(size(xi_C)) beta_C beta_H  mu varthet_1 eta_C alpha_1 beta_CH alpha_2 eta_H xi_H zet_2 zet_1 varthet_2 eta_CH];
[xi_CReg,bintxi_C,rxi_C] = regress(xi_C,Xxi_C);
[RoRegxi_C,bintRoxi_C,rRoxi_C] = regress(Ror,Xxi_C);
% For xi_H
Xxi_H = [ones(size(xi_H)) beta_C beta_H  mu varthet_1 eta_C alpha_1 beta_CH alpha_2 eta_H xi_C zet_2 zet_1 varthet_2 eta_CH];
[xi_HReg,bintxi_H,rxi_H] = regress(xi_H,Xxi_H);
[RoRegxi_H,bintRoxi_H,rRoxi_H] = regress(Ror,Xxi_H);
% For xi_CH
Xxi_CH = [ones(size(xi_CH)) beta_C beta_H  mu varthet_1 eta_C alpha_1 beta_CH alpha_2 eta_H xi_C zet_2 zet_1 varthet_2 eta_CH];
[Psi_KReg,bintxi_CH,rxi_CH] = regress(xi_CH,Xxi_CH);
[RoRegxi_CH,bintRoxi_CH,rRoxi_CH] = regress(Ror,Xxi_CH);
% 
% %We now check for the correlation coefficient, using the Pearson's
% %correlation corr routine in Matlab
rhomu = corr(rmu,rRomu)
rhobeta_C = corr(rbeta_C,rRobeta_C)
rhobeta_H = corr(rbeta_H,rRobeta_H)
rhobeta_CH = corr(rbeta_CH,rRobeta_CH)
rhovarthet_1 = corr(rvarthet_1,rRovarthet_1)
rhovarthet_2 = corr(rvarthet_2,rRovarthet_2)
rhoalpha_1 = corr(ralpha_1,rRoalpha_1)
rhoalpha_2 = corr(ralpha_2,rRoalpha_2)
rhoeta_C = corr(reta_C,rRoeta_C)
rhoeta_H = corr(reta_H,rRoeta_H)
rhoeta_CH = corr(reta_CH,rRoeta_CH)
rhozet_1 = corr(rzet_1,rRozet_1)
rhozet_2 = corr(rzet_2,rRozet_2)
rhoxi_C = corr(rxi_C,rRoxi_C)
rhoxi_H = corr(rxi_H,rRoxi_H)
rhoxi_CH = corr(rxi_CH,rRoxi_CH)
% Now, we want to use the total number of infected persons as the output
% response instead of the reproduction number.

mu1 = D(:,1);
beta_C1 = D(:,2);
beta_H1 = D(:,3);
beta_CH1 = D(:,4);
varthet_11 = D(:,5);
varthet_21 = D(:,6);
alpha_11 = D(:,7);
alpha_21 = D(:,8);
eta_C1 = D(:,9);
eta_H1 = D(:,10);
eta_CH1 = D(:,11);
zet_11 = D(:,12);
zet_21 = D(:,13);
xi_C1 = D(:,14);
xi_H1 = D(:,15);
xi_CH1 = D(:,16);

tspan = [0 10];
%tspan = 0:5/69:5; % see how i break down the spacing for the time i.e the time lenght is broken down into 69 equal spaces between 0 and 5.
yzero = [700000;500;300;500;500;20;100]; % initial conditions.

for i = 1: N % recall what N is.
mu2 = mu1(i);
beta_C2 = beta_C1(i);
beta_H2 = beta_H1(i);
beta_CH2 = beta_CH1(i);
varthet_12 = varthet_11(i);
varthet_22 = varthet_21(i);
alpha_12 = alpha_11(i);
alpha_22 = alpha_21(i);
eta_C2 = eta_C1(i);
eta_H2 = eta_H1(i);
eta_CH2 = eta_CH1(i);
zet_12 = zet_11(i);
zet_22 = zet_21(i);
xi_C2 = xi_C1(i);
xi_H2 = xi_H1(i);
xi_CH2 = xi_CH1(i);
[t,y] = ode45(@SARS_Hepatitis_State,tspan,yzero); %solving the model using each set of parameters from D
Infect(1:length(y(:,2)),i) = y(:,2); %creating a matrix of infected individuals with each set of parameter values
if i == N, break, end
clear t y
end
% Our output response in this case is the total number of infected persons,
% I.
% calculating the output response, which is the total number of infected persons over
%the enrire period of time, for each set of parameter values
for i=1:N
InfectSum(i) = sum(Infect(:,i))
end
% %We now proceed with performing the sensitivity analysis using the PRCC
% %technique and the number of infected persons as response function.
% %We now rank InfectSum
 [InfectSums, InfectSumi] = sort(InfectSum);
 [InfectSumx, InfectSumr] = sort(InfectSumi);
% 
% %InfectSumr is now the rank-transformed vector of InfectSum
Dr;
InfectSumr;
[InfectSumRegmu,bintInfectSummu,rInfectSummu] = regress(InfectSumr',Xmu);
[InfectSumRegbeta_C,bintInfectSumbeta_C,rInfectSumbeta_C] = regress(InfectSumr',Xbeta_C);
[InfectSumRegbeta_H,bintInfectSumbeta_H,rInfectSumbeta_H] = regress(InfectSumr',Xbeta_H);
[InfectSumRegbeta_CH,bintInfectSumbeta_CH,rInfectSumbeta_CH] = regress(InfectSumr',Xbeta_CH);
[InfectSumRegvarthet_1,bintInfectSumvarthet_1,rInfectSumvarthet_1] = regress(InfectSumr',Xvarthet_1);
[InfectSumRegvarthet_2,bintInfectSumvarthet_2,rInfectSumvarthet_2] = regress(InfectSumr',Xvarthet_2);
[InfectSumRegalpha_1,bintInfectSumalpha_1,rInfectSumalpha_1] = regress(InfectSumr',Xalpha_1);
[InfectSumRegalpha_2,bintInfectSumalpha_2,rInfectSumalpha_2] = regress(InfectSumr',Xalpha_2);
[InfectSumRegeta_C,bintInfectSumeta_C,rInfectSumeta_C] = regress(InfectSumr',Xeta_C);
[InfectSumRegeta_H,bintInfectSumeta_H,rInfectSumeta_H] = regress(InfectSumr',Xeta_H);
[InfectSumRegeta_CH,bintInfectSumeta_CH,rInfectSumeta_CH] = regress(InfectSumr',Xeta_CH);
[InfectSumRegzet_1,bintInfectSumzet_1,rInfectSumzet_1] = regress(InfectSumr',Xzet_1);
[InfectSumRegzet_2,bintInfectSumzet_2,rInfectSumzet_2] = regress(InfectSumr',Xzet_2);
[InfectSumRegxi_C,bintInfectSumxi_C,rInfectSumxi_C] = regress(InfectSumr',Xxi_C);
[InfectSumRegxi_H,bintInfectSumxi_H,rInfectSumxi_H] = regress(InfectSumr',Xxi_H);
[InfectSumRegxi_CH,bintInfectSumxi_CH,rInfectSumxi_CH] = regress(InfectSumr',Xxi_CH);
% 
% %We now check for the correlation coefficient, using the Pearson's
% %correlation corr routine when the response function is the number of
% %infected
% 
% 
rho2mu = corr(rmu,rInfectSummu)
rho2beta_C = corr(rbeta_C,rInfectSumbeta_C)
rho2beta_H = corr(rbeta_H,rInfectSumbeta_H)
rho2beta_CH = corr(rbeta_CH,rInfectSumbeta_CH)
rho2varthet_1 = corr(rvarthet_1,rInfectSumvarthet_1)
rho2varthet_2 = corr(rvarthet_2,rInfectSumvarthet_2)
rho2alpha_1 = corr(ralpha_1,rInfectSumalpha_1)
rho2alpha_2 = corr(ralpha_2,rInfectSumalpha_2)
rho2eta_C = corr(reta_C,rInfectSumeta_C)
rho2eta_H = corr(reta_H,rInfectSumeta_H)
rho2eta_CH = corr(reta_CH,rInfectSumeta_CH)
rho2zet_1 = corr(rzet_1,rInfectSumzet_1)
rho2zet_2 = corr(rzet_2,rInfectSumzet_2)
rho2xi_C = corr(rxi_C,rInfectSumxi_C)
rho2xi_H = corr(rxi_H,rInfectSumxi_H)
rho2xi_CH = corr(rxi_CH,rInfectSumxi_CH)

%Bar chat here
figure
C1=[rhomu rhobeta_C rhovarthet_1 rhoalpha_1 rhoxi_C rhoeta_C];
Z1=sort(C1);
bar(C1)
ylabel('PRCC values')
xlabel('Parameters')
set(gca,'XTickLabel',{'\mu','\beta_{C}','\vartheta_{1}','\alpha_{1}','\xi_{C}','\eta_{C}'})
title('Sensitivity analysis with R_{0C} as the response function')
figure
C2=[rho2mu rho2beta_C rho2beta_H rho2beta_CH rho2varthet_1 rho2varthet_2 rho2alpha_1 rho2alpha_2 rho2eta_C rho2eta_H rho2eta_CH rho2zet_1 rho2zet_2 rho2xi_C rho2xi_H rho2xi_CH];
Z2=sort(C2);
bar(C2)
ylabel('PRCC values')
xlabel('Parameters')
set(gca,'XTickLabel',{'\mu','\beta_C','\beta_H','\beta_{CH}','\vartheta_{1}','\vartheta_{2}','\alpha_{1}','\alpha_{2}','\eta_{C}','\eta_{H}','\eta_{CH}','\zeta_{1}','\zeta_{2}','\xi_{C}','\xi_{H}','\xi_{CH}'})
title('Sensitivity analysis with total infected population as the response function')
% toc