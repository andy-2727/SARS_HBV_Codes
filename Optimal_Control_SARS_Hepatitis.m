function Optimal_Control_SARS_Hepatitis
% This function computes the optimal control and the corresponding solution using forward-backward sweep
clc;
clear;
test = -1;
delta = 0.0001; %set tolerance
N = 100; %number of subdivisions
T = 200; %timespan
t=linspace(0,T,N+1);
k2=T/N; %step
k = k2/2;
u1=zeros(1,length(t)); u2=zeros(1,length(t)); u3=zeros(1,length(t)); u4=zeros(1,length(t)); 
S=zeros(1,length(t)); A_C=zeros(1,length(t)); I_C=zeros(1,length(t)); A_H=zeros(1,length(t)); I_H=zeros(1,length(t));
I_CH=zeros(1,length(t)); R=zeros(1,length(t)); 
lambda1=zeros(1,length(t)); lambda2=zeros(1,length(t)); lambda3=zeros(1,length(t)); lambda4=zeros(1,length(t)); lambda5=zeros(1,length(t));
lambda6=zeros(1,length(t)); lambda7=zeros(1,length(t)); 
S(1)=700000; A_C(1)=500; I_C(1)=300; A_H(1)=500; I_H(1)=500; I_CH(1)=20; R(1)=100; 
Lambda = 732500/(76.31*365); mu = 1/(76.31*365); 
beta_C=0.15; beta_H=0.15; beta_CH=0.2; alpha_1 = 1/7; alpha_2 =1/7; eta_C = 0.05; eta_H = 0.05; eta_CH = 0.005; 
varthet_1 = 0.5; varthet_2 =0.5; zet_1 = 1.0; zet_2 = 1.0; xi_C =0.15; xi_H =0.15; xi_CH =0.15;     

R0C = beta_C*(varthet_1*(eta_C + xi_C + mu)+alpha_1)/((alpha_1+mu)*(eta_C + xi_C + mu)) %SARS-CoV-2
R0H = beta_H*(varthet_2*(eta_H + xi_H + mu)+alpha_1)/((alpha_2+mu)*(eta_H + xi_H + mu)) % HBV
RCH = beta_CH/(eta_CH + xi_CH + mu) % co-infection
w1=0.5; w2=0.7; w3=0.9; w4=0.9;  
while (test<0) % while the tolerance is reached, repeat
oldu1 = u1;
oldu2 = u2;
oldu3 = u3;
oldu4 = u4;
oldS = S;
oldA_C = A_C;
oldI_C = I_C;
oldA_H = A_H;
oldI_H = I_H;
oldI_CH = I_CH;
oldR = R;
oldlambda1 = lambda1;
oldlambda2 = lambda2;
oldlambda3 = lambda3;
oldlambda4 = lambda4;
oldlambda5 = lambda5;
oldlambda6 = lambda6;
oldlambda7 = lambda7;
for i=1:N %loop that solve the forward differential equation of the state system
auxS1= Lambda - (1-u1(i))*beta_C*(varthet_1*A_C(i) + I_C(i))/(S(i) + A_C(i) + I_C(i) + A_H(i) + I_H(i) + I_CH(i) + R(i))*S(i) - (1-u2(i))*beta_H*(varthet_2*A_H(i) + I_H(i))/(S(i) + A_C(i) + I_C(i) + A_H(i) + I_H(i) + I_CH(i) + R(i))*S(i) - (1-u3(i))*(beta_CH*I_CH(i))/(S(i) + A_C(i) + I_C(i) + A_H(i) + I_H(i) + I_CH(i) + R(i))*S(i) - mu*S(i); 
auxA_C1 = (1-u1(i))*beta_C*(varthet_1*A_C(i) + I_C(i))/(S(i) + A_C(i) + I_C(i) + A_H(i) + I_H(i) + I_CH(i) + R(i))*S(i) - (alpha_1 + mu)*A_C(i); 
auxI_C1 = alpha_1*A_C(i) - (xi_C + eta_C + mu)*I_C(i) - (1-u4(i))*zet_1*beta_H*(varthet_2*A_H(i) + I_H(i))/(S(i) + A_C(i) + I_C(i) + A_H(i) + I_H(i) + I_CH(i) + R(i))*I_C(i);
auxA_H1 = (1-u2(i))*beta_H*(varthet_2*A_H(i) + I_H(i))/(S(i) + A_C(i) + I_C(i) + A_H(i) + I_H(i) + I_CH(i) + R(i))*S(i) - (alpha_2 + mu)*A_H(i); 
auxI_H1 = alpha_2*A_H(i) - (xi_H + eta_H + mu)*I_H(i) - (1-u4(i))*zet_2*beta_C*(varthet_1*A_C(i) + I_C(i))/(S(i) + A_C(i) + I_C(i) + A_H(i) + I_H(i) + I_CH(i) + R(i))*I_H(i); 
auxI_CH1 = (1-u3(i))*(beta_CH*I_CH(i))/(S(i) + A_C(i) + I_C(i) + A_H(i) + I_H(i) + I_CH(i) + R(i))*S(i) + (1-u4(i))*zet_1*beta_H*(varthet_2*A_H(i) + I_H(i))/(S(i) + A_C(i) + I_C(i) + A_H(i) + I_H(i) + I_CH(i) + R(i))*I_C(i) + (1-u4(i))*zet_2*beta_C*(varthet_1*A_C(i) + I_C(i))/(S(i) + A_C(i) + I_C(i) + A_H(i) + I_H(i) + I_CH(i) + R(i))*I_H(i) - (xi_CH + eta_CH + mu)*I_CH(i); 
auxR1 = xi_C*I_C(i) + xi_H*I_H(i) + xi_CH*I_CH(i) - mu*R(i);
S(i+1) = S(i) + k*auxS1;
A_C(i+1) = A_C(i) + k*auxA_C1;  
I_C(i+1) = I_C(i) + k*auxI_C1;
A_H(i+1) = A_H(i) + k*auxA_H1;
I_H(i+1) = I_H(i) + k*auxI_H1;
I_CH(i+1) = I_CH(i) + k*auxI_CH1; 
R(i+1)    =  R(i) + k*auxR1;
end
for i=1:N %loop that solves the backward differential equation of the adjoint system
j = N + 2 -i;
auxlambda11 = lambda1(j)*(mu - (beta_C*(I_C(j) + A_C(j)*varthet_1)*(u1(j) - 1))/(A_C(j) + A_H(j) + I_C(j) + I_CH(j) + I_H(j) + R(j) + S(j)) - (beta_H*(I_H(j) + A_H(j)*varthet_2)*(u2(j) - 1))/(A_C(j) + A_H(j) + I_C(j) + I_CH(j) + I_H(j) + R(j) + S(j)) - (I_CH(j)*beta_CH*(u3(j) - 1))/(A_C(j) + A_H(j) + I_C(j) + I_CH(j) + I_H(j) + R(j) + S(j)) + (S(j)*beta_C*(I_C(j) + A_C(j)*varthet_1)*(u1(j) - 1))/(A_C(j) + A_H(j) + I_C(j) + I_CH(j) + I_H(j) + R(j) + S(j))^2 + (S(j)*beta_H*(I_H(j) + A_H(j)*varthet_2)*(u2(j) - 1))/(A_C(j) + A_H(j) + I_C(j) + I_CH(j) + I_H(j) + R(j) + S(j))^2 + (I_CH(j)*S(j)*beta_CH*(u3(j) - 1))/(A_C(j) + A_H(j) + I_C(j) + I_CH(j) + I_H(j) + R(j) + S(j))^2) - lambda6(j)*((I_CH(j)*S(j)*beta_CH*(u3(j) - 1))/(A_C(j) + A_H(j) + I_C(j) + I_CH(j) + I_H(j) + R(j) + S(j))^2 - (I_CH(j)*beta_CH*(u3(j) - 1))/(A_C(j) + A_H(j) + I_C(j) + I_CH(j) + I_H(j) + R(j) + S(j)) + (I_H(j)*beta_C*zet_2*(I_C(j) + A_C(j)*varthet_1)*(u4(j) - 1))/(A_C(j) + A_H(j) + I_C(j) + I_CH(j) + I_H(j) + R(j) + S(j))^2 + (I_C(j)*beta_H*zet_1*(I_H(j) + A_H(j)*varthet_2)*(u4(j) - 1))/(A_C(j) + A_H(j) + I_C(j) + I_CH(j) + I_H(j) + R(j) + S(j))^2) + lambda2(j)*((beta_C*(I_C(j) + A_C(j)*varthet_1)*(u1(j) - 1))/(A_C(j) + A_H(j) + I_C(j) + I_CH(j) + I_H(j) + R(j) + S(j)) - (S(j)*beta_C*(I_C(j) + A_C(j)*varthet_1)*(u1(j) - 1))/(A_C(j) + A_H(j) + I_C(j) + I_CH(j) + I_H(j) + R(j) + S(j))^2) + lambda4(j)*((beta_H*(I_H(j) + A_H(j)*varthet_2)*(u2(j) - 1))/(A_C(j) + A_H(j) + I_C(j) + I_CH(j) + I_H(j) + R(j) + S(j)) - (S(j)*beta_H*(I_H(j) + A_H(j)*varthet_2)*(u2(j) - 1))/(A_C(j) + A_H(j) + I_C(j) + I_CH(j) + I_H(j) + R(j) + S(j))^2) + (I_H(j)*beta_C*lambda5(j)*zet_2*(I_C(j) + A_C(j)*varthet_1)*(u4(j) - 1))/(A_C(j) + A_H(j) + I_C(j) + I_CH(j) + I_H(j) + R(j) + S(j))^2 + (I_C(j)*beta_H*lambda3(j)*zet_1*(I_H(j) + A_H(j)*varthet_2)*(u4(j) - 1))/(A_C(j) + A_H(j) + I_C(j) + I_CH(j) + I_H(j) + R(j) + S(j))^2;
auxlambda21 = lambda1(j)*((S(j)*beta_C*(I_C(j) + A_C(j)*varthet_1)*(u1(j) - 1))/(A_C(j) + A_H(j) + I_C(j) + I_CH(j) + I_H(j) + R(j) + S(j))^2 - (S(j)*beta_C*varthet_1*(u1(j) - 1))/(A_C(j) + A_H(j) + I_C(j) + I_CH(j) + I_H(j) + R(j) + S(j)) + (S(j)*beta_H*(I_H(j) + A_H(j)*varthet_2)*(u2(j) - 1))/(A_C(j) + A_H(j) + I_C(j) + I_CH(j) + I_H(j) + R(j) + S(j))^2 + (I_CH(j)*S(j)*beta_CH*(u3(j) - 1))/(A_C(j) + A_H(j) + I_C(j) + I_CH(j) + I_H(j) + R(j) + S(j))^2) - lambda3(j)*(alpha_1 - (I_C(j)*beta_H*zet_1*(I_H(j) + A_H(j)*varthet_2)*(u4(j) - 1))/(A_C(j) + A_H(j) + I_C(j) + I_CH(j) + I_H(j) + R(j) + S(j))^2) - lambda6(j)*((I_CH(j)*S(j)*beta_CH*(u3(j) - 1))/(A_C(j) + A_H(j) + I_C(j) + I_CH(j) + I_H(j) + R(j) + S(j))^2 - (I_H(j)*beta_C*varthet_1*zet_2*(u4(j) - 1))/(A_C(j) + A_H(j) + I_C(j) + I_CH(j) + I_H(j) + R(j) + S(j)) + (I_H(j)*beta_C*zet_2*(I_C(j) + A_C(j)*varthet_1)*(u4(j) - 1))/(A_C(j) + A_H(j) + I_C(j) + I_CH(j) + I_H(j) + R(j) + S(j))^2 + (I_C(j)*beta_H*zet_1*(I_H(j) + A_H(j)*varthet_2)*(u4(j) - 1))/(A_C(j) + A_H(j) + I_C(j) + I_CH(j) + I_H(j) + R(j) + S(j))^2) + lambda2(j)*(alpha_1 + mu + (S(j)*beta_C*varthet_1*(u1(j) - 1))/(A_C(j) + A_H(j) + I_C(j) + I_CH(j) + I_H(j) + R(j) + S(j)) - (S(j)*beta_C*(I_C(j) + A_C(j)*varthet_1)*(u1(j) - 1))/(A_C(j) + A_H(j) + I_C(j) + I_CH(j) + I_H(j) + R(j) + S(j))^2) - lambda5(j)*((I_H(j)*beta_C*varthet_1*zet_2*(u4(j) - 1))/(A_C(j) + A_H(j) + I_C(j) + I_CH(j) + I_H(j) + R(j) + S(j)) - (I_H(j)*beta_C*zet_2*(I_C(j) + A_C(j)*varthet_1)*(u4(j) - 1))/(A_C(j) + A_H(j) + I_C(j) + I_CH(j) + I_H(j) + R(j) + S(j))^2) - (S(j)*beta_H*lambda4(j)*(I_H(j) + A_H(j)*varthet_2)*(u2(j) - 1))/(A_C(j) + A_H(j) + I_C(j) + I_CH(j) + I_H(j) + R(j) + S(j))^2 - 1;
auxlambda31 = lambda2(j)*((S(j)*beta_C*(u1(j) - 1))/(A_C(j) + A_H(j) + I_C(j) + I_CH(j) + I_H(j) + R(j) + S(j)) - (S(j)*beta_C*(I_C(j) + A_C(j)*varthet_1)*(u1(j) - 1))/(A_C(j) + A_H(j) + I_C(j) + I_CH(j) + I_H(j) + R(j) + S(j))^2) - lambda5(j)*((I_H(j)*beta_C*zet_2*(u4(j) - 1))/(A_C(j) + A_H(j) + I_C(j) + I_CH(j) + I_H(j) + R(j) + S(j)) - (I_H(j)*beta_C*zet_2*(I_C(j) + A_C(j)*varthet_1)*(u4(j) - 1))/(A_C(j) + A_H(j) + I_C(j) + I_CH(j) + I_H(j) + R(j) + S(j))^2) - lambda7(j)*xi_C - lambda6(j)*((I_CH(j)*S(j)*beta_CH*(u3(j) - 1))/(A_C(j) + A_H(j) + I_C(j) + I_CH(j) + I_H(j) + R(j) + S(j))^2 - (beta_H*zet_1*(I_H(j) + A_H(j)*varthet_2)*(u4(j) - 1))/(A_C(j) + A_H(j) + I_C(j) + I_CH(j) + I_H(j) + R(j) + S(j)) - (I_H(j)*beta_C*zet_2*(u4(j) - 1))/(A_C(j) + A_H(j) + I_C(j) + I_CH(j) + I_H(j) + R(j) + S(j)) + (I_H(j)*beta_C*zet_2*(I_C(j) + A_C(j)*varthet_1)*(u4(j) - 1))/(A_C(j) + A_H(j) + I_C(j) + I_CH(j) + I_H(j) + R(j) + S(j))^2 + (I_C(j)*beta_H*zet_1*(I_H(j) + A_H(j)*varthet_2)*(u4(j) - 1))/(A_C(j) + A_H(j) + I_C(j) + I_CH(j) + I_H(j) + R(j) + S(j))^2) + lambda1(j)*((S(j)*beta_C*(I_C(j) + A_C(j)*varthet_1)*(u1(j) - 1))/(A_C(j) + A_H(j) + I_C(j) + I_CH(j) + I_H(j) + R(j) + S(j))^2 - (S(j)*beta_C*(u1(j) - 1))/(A_C(j) + A_H(j) + I_C(j) + I_CH(j) + I_H(j) + R(j) + S(j)) + (S(j)*beta_H*(I_H(j) + A_H(j)*varthet_2)*(u2(j) - 1))/(A_C(j) + A_H(j) + I_C(j) + I_CH(j) + I_H(j) + R(j) + S(j))^2 + (I_CH(j)*S(j)*beta_CH*(u3(j) - 1))/(A_C(j) + A_H(j) + I_C(j) + I_CH(j) + I_H(j) + R(j) + S(j))^2) + lambda3(j)*(eta_C + mu + xi_C - (beta_H*zet_1*(I_H(j) + A_H(j)*varthet_2)*(u4(j) - 1))/(A_C(j) + A_H(j) + I_C(j) + I_CH(j) + I_H(j) + R(j) + S(j)) + (I_C(j)*beta_H*zet_1*(I_H(j) + A_H(j)*varthet_2)*(u4(j) - 1))/(A_C(j) + A_H(j) + I_C(j) + I_CH(j) + I_H(j) + R(j) + S(j))^2) - (S(j)*beta_H*lambda4(j)*(I_H(j) + A_H(j)*varthet_2)*(u2(j) - 1))/(A_C(j) + A_H(j) + I_C(j) + I_CH(j) + I_H(j) + R(j) + S(j))^2 - 1;
auxlambda41 = lambda1(j)*((S(j)*beta_C*(I_C(j) + A_C(j)*varthet_1)*(u1(j) - 1))/(A_C(j) + A_H(j) + I_C(j) + I_CH(j) + I_H(j) + R(j) + S(j))^2 - (S(j)*beta_H*varthet_2*(u2(j) - 1))/(A_C(j) + A_H(j) + I_C(j) + I_CH(j) + I_H(j) + R(j) + S(j)) + (S(j)*beta_H*(I_H(j) + A_H(j)*varthet_2)*(u2(j) - 1))/(A_C(j) + A_H(j) + I_C(j) + I_CH(j) + I_H(j) + R(j) + S(j))^2 + (I_CH(j)*S(j)*beta_CH*(u3(j) - 1))/(A_C(j) + A_H(j) + I_C(j) + I_CH(j) + I_H(j) + R(j) + S(j))^2) - lambda5(j)*(alpha_2 - (I_H(j)*beta_C*zet_2*(I_C(j) + A_C(j)*varthet_1)*(u4(j) - 1))/(A_C(j) + A_H(j) + I_C(j) + I_CH(j) + I_H(j) + R(j) + S(j))^2) - lambda6(j)*((I_CH(j)*S(j)*beta_CH*(u3(j) - 1))/(A_C(j) + A_H(j) + I_C(j) + I_CH(j) + I_H(j) + R(j) + S(j))^2 - (I_C(j)*beta_H*varthet_2*zet_1*(u4(j) - 1))/(A_C(j) + A_H(j) + I_C(j) + I_CH(j) + I_H(j) + R(j) + S(j)) + (I_H(j)*beta_C*zet_2*(I_C(j) + A_C(j)*varthet_1)*(u4(j) - 1))/(A_C(j) + A_H(j) + I_C(j) + I_CH(j) + I_H(j) + R(j) + S(j))^2 + (I_C(j)*beta_H*zet_1*(I_H(j) + A_H(j)*varthet_2)*(u4(j) - 1))/(A_C(j) + A_H(j) + I_C(j) + I_CH(j) + I_H(j) + R(j) + S(j))^2) + lambda4(j)*(alpha_2 + mu + (S(j)*beta_H*varthet_2*(u2(j) - 1))/(A_C(j) + A_H(j) + I_C(j) + I_CH(j) + I_H(j) + R(j) + S(j)) - (S(j)*beta_H*(I_H(j) + A_H(j)*varthet_2)*(u2(j) - 1))/(A_C(j) + A_H(j) + I_C(j) + I_CH(j) + I_H(j) + R(j) + S(j))^2) - lambda3(j)*((I_C(j)*beta_H*varthet_2*zet_1*(u4(j) - 1))/(A_C(j) + A_H(j) + I_C(j) + I_CH(j) + I_H(j) + R(j) + S(j)) - (I_C(j)*beta_H*zet_1*(I_H(j) + A_H(j)*varthet_2)*(u4(j) - 1))/(A_C(j) + A_H(j) + I_C(j) + I_CH(j) + I_H(j) + R(j) + S(j))^2) - (S(j)*beta_C*lambda2(j)*(I_C(j) + A_C(j)*varthet_1)*(u1(j) - 1))/(A_C(j) + A_H(j) + I_C(j) + I_CH(j) + I_H(j) + R(j) + S(j))^2 - 1;
auxlambda51 = lambda4(j)*((S(j)*beta_H*(u2(j) - 1))/(A_C(j) + A_H(j) + I_C(j) + I_CH(j) + I_H(j) + R(j) + S(j)) - (S(j)*beta_H*(I_H(j) + A_H(j)*varthet_2)*(u2(j) - 1))/(A_C(j) + A_H(j) + I_C(j) + I_CH(j) + I_H(j) + R(j) + S(j))^2) - lambda3(j)*((I_C(j)*beta_H*zet_1*(u4(j) - 1))/(A_C(j) + A_H(j) + I_C(j) + I_CH(j) + I_H(j) + R(j) + S(j)) - (I_C(j)*beta_H*zet_1*(I_H(j) + A_H(j)*varthet_2)*(u4(j) - 1))/(A_C(j) + A_H(j) + I_C(j) + I_CH(j) + I_H(j) + R(j) + S(j))^2) - lambda7(j)*xi_H - lambda6(j)*((I_CH(j)*S(j)*beta_CH*(u3(j) - 1))/(A_C(j) + A_H(j) + I_C(j) + I_CH(j) + I_H(j) + R(j) + S(j))^2 - (beta_C*zet_2*(I_C(j) + A_C(j)*varthet_1)*(u4(j) - 1))/(A_C(j) + A_H(j) + I_C(j) + I_CH(j) + I_H(j) + R(j) + S(j)) - (I_C(j)*beta_H*zet_1*(u4(j) - 1))/(A_C(j) + A_H(j) + I_C(j) + I_CH(j) + I_H(j) + R(j) + S(j)) + (I_H(j)*beta_C*zet_2*(I_C(j) + A_C(j)*varthet_1)*(u4(j) - 1))/(A_C(j) + A_H(j) + I_C(j) + I_CH(j) + I_H(j) + R(j) + S(j))^2 + (I_C(j)*beta_H*zet_1*(I_H(j) + A_H(j)*varthet_2)*(u4(j) - 1))/(A_C(j) + A_H(j) + I_C(j) + I_CH(j) + I_H(j) + R(j) + S(j))^2) + lambda1(j)*((S(j)*beta_C*(I_C(j) + A_C(j)*varthet_1)*(u1(j) - 1))/(A_C(j) + A_H(j) + I_C(j) + I_CH(j) + I_H(j) + R(j) + S(j))^2 - (S(j)*beta_H*(u2(j) - 1))/(A_C(j) + A_H(j) + I_C(j) + I_CH(j) + I_H(j) + R(j) + S(j)) + (S(j)*beta_H*(I_H(j) + A_H(j)*varthet_2)*(u2(j) - 1))/(A_C(j) + A_H(j) + I_C(j) + I_CH(j) + I_H(j) + R(j) + S(j))^2 + (I_CH(j)*S(j)*beta_CH*(u3(j) - 1))/(A_C(j) + A_H(j) + I_C(j) + I_CH(j) + I_H(j) + R(j) + S(j))^2) + lambda5(j)*(eta_H + mu + xi_H - (beta_C*zet_2*(I_C(j) + A_C(j)*varthet_1)*(u4(j) - 1))/(A_C(j) + A_H(j) + I_C(j) + I_CH(j) + I_H(j) + R(j) + S(j)) + (I_H(j)*beta_C*zet_2*(I_C(j) + A_C(j)*varthet_1)*(u4(j) - 1))/(A_C(j) + A_H(j) + I_C(j) + I_CH(j) + I_H(j) + R(j) + S(j))^2) - (S(j)*beta_C*lambda2(j)*(I_C(j) + A_C(j)*varthet_1)*(u1(j) - 1))/(A_C(j) + A_H(j) + I_C(j) + I_CH(j) + I_H(j) + R(j) + S(j))^2 - 1;
auxlambda61 = lambda1(j)*((S(j)*beta_C*(I_C(j) + A_C(j)*varthet_1)*(u1(j) - 1))/(A_C(j) + A_H(j) + I_C(j) + I_CH(j) + I_H(j) + R(j) + S(j))^2 - (S(j)*beta_CH*(u3(j) - 1))/(A_C(j) + A_H(j) + I_C(j) + I_CH(j) + I_H(j) + R(j) + S(j)) + (S(j)*beta_H*(I_H(j) + A_H(j)*varthet_2)*(u2(j) - 1))/(A_C(j) + A_H(j) + I_C(j) + I_CH(j) + I_H(j) + R(j) + S(j))^2 + (I_CH(j)*S(j)*beta_CH*(u3(j) - 1))/(A_C(j) + A_H(j) + I_C(j) + I_CH(j) + I_H(j) + R(j) + S(j))^2) - lambda7(j)*xi_CH + lambda6(j)*(eta_CH +  mu + xi_CH + (S(j)*beta_CH*(u3(j) - 1))/(A_C(j) + A_H(j) + I_C(j) + I_CH(j) + I_H(j) + R(j) + S(j)) - (I_CH(j)*S(j)*beta_CH*(u3(j) - 1))/(A_C(j) + A_H(j) + I_C(j) + I_CH(j) + I_H(j) + R(j) + S(j))^2 - (I_H(j)*beta_C*zet_2*(I_C(j) + A_C(j)*varthet_1)*(u4(j) - 1))/(A_C(j) + A_H(j) + I_C(j) + I_CH(j) + I_H(j) + R(j) + S(j))^2 - (I_C(j)*beta_H*zet_1*(I_H(j) + A_H(j)*varthet_2)*(u4(j) - 1))/(A_C(j) + A_H(j) + I_C(j) + I_CH(j) + I_H(j) + R(j) + S(j))^2) - (S(j)*beta_C*lambda2(j)*(I_C(j) + A_C(j)*varthet_1)*(u1(j) - 1))/(A_C(j) + A_H(j) + I_C(j) + I_CH(j) + I_H(j) + R(j) + S(j))^2 - (S(j)*beta_H*lambda4(j)*(I_H(j) + A_H(j)*varthet_2)*(u2(j) - 1))/(A_C(j) + A_H(j) + I_C(j) + I_CH(j) + I_H(j) + R(j) + S(j))^2 + (I_H(j)*beta_C*lambda5(j)*zet_2*(I_C(j) + A_C(j)*varthet_1)*(u4(j) - 1))/(A_C(j) + A_H(j) + I_C(j) + I_CH(j) + I_H(j) + R(j) + S(j))^2 + (I_C(j)*beta_H*lambda3(j)*zet_1*(I_H(j) + A_H(j)*varthet_2)*(u4(j) - 1))/(A_C(j) + A_H(j) + I_C(j) + I_CH(j) + I_H(j) + R(j) + S(j))^2 - 1;
auxlambda71 = lambda7(j)*mu - lambda6(j)*((I_CH(j)*S(j)*beta_CH*(u3(j) - 1))/(A_C(j) + A_H(j) + I_C(j) + I_CH(j) + I_H(j) + R(j) + S(j))^2 + (I_H(j)*beta_C*zet_2*(I_C(j) + A_C(j)*varthet_1)*(u4(j) - 1))/(A_C(j) + A_H(j) + I_C(j) + I_CH(j) + I_H(j) + R(j) + S(j))^2 + (I_C(j)*beta_H*zet_1*(I_H(j) + A_H(j)*varthet_2)*(u4(j) - 1))/(A_C(j) + A_H(j) + I_C(j) + I_CH(j) + I_H(j) + R(j) + S(j))^2) + lambda1(j)*((S(j)*beta_C*(I_C(j) + A_C(j)*varthet_1)*(u1(j) - 1))/(A_C(j) + A_H(j) + I_C(j) + I_CH(j) + I_H(j) + R(j) + S(j))^2 + (S(j)*beta_H*(I_H(j) + A_H(j)*varthet_2)*(u2(j) - 1))/(A_C(j) + A_H(j) + I_C(j) + I_CH(j) + I_H(j) + R(j) + S(j))^2 + (I_CH(j)*S(j)*beta_CH*(u3(j) - 1))/(A_C(j) + A_H(j) + I_C(j) + I_CH(j) + I_H(j) + R(j) + S(j))^2) - (S(j)*beta_C*lambda2(j)*(I_C(j) + A_C(j)*varthet_1)*(u1(j) - 1))/(A_C(j) + A_H(j) + I_C(j) + I_CH(j) + I_H(j) + R(j) + S(j))^2 - (S(j)*beta_H*lambda4(j)*(I_H(j) + A_H(j)*varthet_2)*(u2(j) - 1))/(A_C(j) + A_H(j) + I_C(j) + I_CH(j) + I_H(j) + R(j) + S(j))^2 + (I_H(j)*beta_C*lambda5(j)*zet_2*(I_C(j) + A_C(j)*varthet_1)*(u4(j) - 1))/(A_C(j) + A_H(j) + I_C(j) + I_CH(j) + I_H(j) + R(j) + S(j))^2 + (I_C(j)*beta_H*lambda3(j)*zet_1*(I_H(j) + A_H(j)*varthet_2)*(u4(j) - 1))/(A_C(j) + A_H(j) + I_C(j) + I_CH(j) + I_H(j) + R(j) + S(j))^2; 
lambda1(j-1) = lambda1(j) - k*auxlambda11;
lambda2(j-1) = lambda2(j) - k*auxlambda21;
lambda3(j-1) = lambda3(j) - k*auxlambda31; 
lambda4(j-1) = lambda4(j) - k*auxlambda41;
lambda5(j-1) = lambda5(j) - k*auxlambda51;
lambda6(j-1) = lambda6(j) - k*auxlambda61;
lambda7(j-1) = lambda7(j) - k*auxlambda71;
end
v1Aux= (1/w1)*(beta_C.*(varthet_1.*A_C+I_C)./(S + A_C + I_C + A_H + I_H + I_CH + R)).*(S.*(lambda2-lambda1));
v2Aux= (1/w2)*(beta_H.*(varthet_2.*A_H+I_H)./(S + A_C + I_C + A_H + I_H + I_CH + R)).*(S.*(lambda4-lambda1));
v3Aux= (1/w3)*(beta_CH.*I_CH/(S + A_C + I_C + A_H + I_H + I_CH + R)).*(S.*(lambda6-lambda1));
v4Aux= (1/w4)*((beta_C.*zet_2.*(varthet_1.*A_C+I_C)./(S + A_C + I_C + A_H + I_H + I_CH + R)).*(I_H.*(lambda6-lambda5)) ...
       + (beta_H.*zet_1.*(varthet_2.*A_H+I_H)./(S + A_C + I_C + A_H + I_H + I_CH + R)).*(I_C.*(lambda6-lambda3)));
auxu1 = min(1,max(0,v1Aux));
auxu2 = min(1,max(0,v2Aux));
auxu3 = min(1,max(0,v3Aux));
auxu4 = min(1,max(0,v4Aux));
u1 = 0.5*(auxu1 + oldu1);
u2 = 0.0*(auxu2 + oldu2);
u3 = 0.0*(auxu3 + oldu3);
u4 = 0.0*(auxu4 + oldu4);
temp1 = delta*sum(abs(u1)) - sum(abs(oldu1 - u1));
temp2 = delta*sum(abs(u2)) - sum(abs(oldu2 - u2));
temp3 = delta*sum(abs(u3)) - sum(abs(oldu3 - u3));
temp4 = delta*sum(abs(u4)) - sum(abs(oldu4 - u4));
temp5 = delta*sum(abs(S)) - sum(abs(oldS - S));
temp6 = delta*sum(abs(A_C)) - sum(abs(oldA_C - A_C));
temp7 = delta*sum(abs(I_C)) - sum(abs(oldI_C - I_C));
temp8 = delta*sum(abs(A_H)) - sum(abs(oldA_H - A_H));
temp9 = delta*sum(abs(I_H)) - sum(abs(oldI_H - I_H));
temp10 = delta*sum(abs(I_CH)) - sum(abs(oldI_CH - I_CH));
temp11 = delta*sum(abs(R)) - sum(abs(oldR - R));
temp12 = delta*sum(abs(lambda1)) - sum(abs(oldlambda1 - lambda1));
temp13 = delta*sum(abs(lambda2)) - sum(abs(oldlambda2 - lambda2));
temp14 = delta*sum(abs(lambda3)) - sum(abs(oldlambda3 - lambda3));
temp15 = delta*sum(abs(lambda4)) - sum(abs(oldlambda4 - lambda4));
temp16 = delta*sum(abs(lambda5)) - sum(abs(oldlambda5 - lambda5));
temp17 = delta*sum(abs(lambda6)) - sum(abs(oldlambda6 - lambda6));
temp18 = delta*sum(abs(lambda7)) - sum(abs(oldlambda7 - lambda7));
test1 = min(temp10, min(temp11, min(temp12, min(temp13, min(temp14, min(temp15, min(temp16, min(temp17, min(temp18)))))))));
test = min(temp1, min(temp2, min(temp3, min(temp4, min(temp5, min(temp6, min(temp7, min(temp8, min(temp9, test1)))))))));
end
plot(t,u1,'b','linewidth',5)
xlabel('Time (days)'),ylabel('Control profile')
legend('u_1')
%plot(t,I_CH,'r','linewidth',5)
%xlabel('Time (days)'),ylabel('I_{CH}')
hold on
