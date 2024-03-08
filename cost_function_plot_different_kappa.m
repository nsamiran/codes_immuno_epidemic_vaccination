% clc
% clear all
% close all 
format long
set(0,'DefaultAxesFontSize',15);

N=10000000;
dt=0.05;%0.01;
%t1=4*30;
TTT=2000;
t=0:dt:TTT;
L=0.75;d=0.001;%interesting 0.003;
T=80.73;V0=500;%d=0.433;T=24.73;
% V=L*N./(1+exp(-d*(t-T)));
%V1=(L*N*d*exp(-d*(t-T)))./((1+exp(-d*(t-T))).^2);
V=L*(N-(N-V0)*exp(-d*t));
V1=L*d*(N-V0)*exp(-d*t);
% phi=0.9411*exp(-((t-117.8)/92.44).^2);
phi=0.01152*t.^(1.023).*exp(-5.01*10^(-6)*t.^(2.412));
psi=1.034*exp(-((t+206.6)/1133).^2);
viral_load=(1.829*10^5)*exp(-((t-3.136)/1.294).^2);
c=0.5*10^(-5);%1.02*10^(-5);
beta=c*viral_load;
alpha=0.8;
kappa=1;
b=0.2;
eps=0.7;


for i=1:length(t)
  
    MM(i)=(dt/2)*(phi(i)*V1(1)+phi(1)*V1(i)+2*sum(phi(i-1:-1:2).*V1(2:1:i-1)));
     MM1(i)=(alpha/N)*MM(i);
          MM2(i)=((1-alpha)/N)*MM(i);

end

% % % % % a1=34.55447;b1=0.60847;aa1=226.40545;bb1=0.17171;a2=35.00855;b2=0.60511;aa2=186.11379;bb2=0.20636;
a1=32.17136;b1=0.2206;aa1=65.40545;bb1=0.210; a2=36.02855;b2=0.57511;aa2=40.11379;bb2=0.27636;
A1= 0.85/((b1^a1)* gamma(a1));AA1= 0.15/((bb1^aa1)* gamma(aa1));A2= 0.94/((b2^a2)* gamma(a2));AA2= 0.06/((bb2^aa2)* gamma(aa2));p1=.9975;p2=1-p1;

%%%%%%%%%%%%%%%%%%%% Probability density model %%%%%%%%%%%%%%%%%%%%

t(1)=0;
S(1)=N-1;
I(1)=1;%data7d(1);%249361;
P(1)=0;
R(1)=0;
D(1)=0;
J1(1)=1;
J2(1)=0;
J(1)=J1(1)+J2(1);
CC=0;
ngap=3:1:30;
for kk=1:length(ngap)
t(1)=0;
t1=ngap(kk)*30; gap=floor(TTT/t1)+1;ngap1(kk)=floor(TTT/t1)+1;
T1=t1/dt;
CJ(1)=J(1);

for m=1:gap-1
for j=1+T1*(m-1):T1*m
      t(j+1)=t(j)+dt;
    CC=CC+1
  %  t(j+1)=t(j)+dt;
    h=dt;
    
    %%%%%%%%%%%beta integration%%%%%%%%%%

binte1(j)=(dt/2)*(beta(1)*J(j)+ beta(j)*J(1));
    for k=2 :j-1
        Y1(k)=dt*beta(k)*J(j-k+1);
    binte1(j)=binte1(j) + Y1(k);
    end
    
    BETA(j)=binte1(j);
    
    
    %%%%%%%%%%%recovery integration%%%%%%%%%%
    inte1(j)=(dt/2)*((exp(-t(j)/b1))*((t(j))^(a1-1))*J(1)+0);
    for k=2 :j-1
        Y1(k)=dt*(J(k))*(exp(-t(j-k+1)/b1))*(t(j-k+1))^(a1-1);
    inte1(j)=inte1(j) + Y1(k);
    end
    
   inte11(j)=(dt/2)*((exp(-t(j)/bb1))*((t(j))^(aa1-1))*J(1)+0);
    for k=2 :j-1
        Y11(k)=dt*(J(k))*(exp(-t(j-k+1)/bb1))*(t(j-k+1))^(aa1-1);
    inte11(j)=inte11(j) + Y11(k);
    end
    
    %%%%%%%%%%%dead integration%%%%%%%%%%
    inte2(j)=(dt/2)*((exp(-t(j)/b2))*(t(j))^(a2-1)*J(1)+0);
    for k=2 :j-1
      Y2(k)=dt*(J(k))*(exp(-t(j-k+1)/b2))*(t(j-k+1))^(a2-1);
    inte2(j)=inte2(j) + Y2(k);
    end
    
    inte22(j)=(dt/2)*((exp(-t(j)/bb2))*(t(j))^(aa2-1)*J(1)+0);
    for k=2 :j-1
      Y22(k)=dt*(J(k))*(exp(-t(j-k+1)/bb2))*(t(j-k+1))^(aa2-1);
    inte22(j)=inte22(j) + Y22(k);
    end
    
    
    AA(j)=p1*A1*inte1(j)+p1*AA1*inte11(j)+p2*A2*inte2(j)+p2*AA2*inte22(j);
    
    %%%%%%%%%%%level of immunity integration%%%%%%%%%%


    
    
      Ms(j)=(dt/2)*(psi(j)*(p1*A1*inte1(1)+p1*AA1*inte11(1))+psi(1)*(p1*A1*inte1(j)+p1*AA1*inte11(j)));%+2*(sum(psi(i-1:-1:2).*recovery_data(2:1:i-1))));
          
           for k=2:j-1
    
                 Ms(j)=Ms(j)+dt*psi(j-k+1)*(p1*A1*inte1(k)+p1*AA1*inte11(k));
           end
            M1(j)=(1-b)*(Ms(j)/N);
           M2(j)=eps*b*(Ms(j)/N);
           MV1(j)=0;
           MV2(j)=0;
           for n=1:m
           M1(j)=M1(j)+MM1(j-T1*(n-1));
           M2(j)=M2(j)+MM2(j-T1*(n-1));
            MV1(j)=MV1(j)+MM1(j-T1*(n-1));
           MV2(j)=MV2(j)+MM2(j-T1*(n-1));
           end
           M(j)=M1(j)+M2(j);
            MV(j)=MV1(j)+MV2(j);
           RR(j)=(dt/2)*(R(1)+R(j)+2*sum(R(2:j-1)));
           AAA(j)=dt*(p1*A1*inte1(j)+p1*AA1*inte11(j));
             RRR(j)=(dt/2)*(AAA(1)+AAA(j)+2*sum(AAA(2:j-1)));
           JJ2(j)=(dt/2)*(J2(1)+J2(j)+2*sum(J2(2:j-1)));
   % Mf(j)=MM1(j)+(Ms(j)/N);  
%      S(j+1)=S(j)+dt*(-BETA(j)*(1-Mf(j)-((I(j)+D(j))/N)));
    S(j+1)=N-(I(j)+D(j)+P(j)+M1(j)*N);
    P(j+1)=b*RRR(j)-JJ2(j)-M2(j)*N;
    J1(j+1)=(BETA(j)*S(j)/N);
    J2(j+1)=(kappa*BETA(j)*P(j)/N);
    J(j+1)=J1(j+1)+J2(j+1);
    I(j+1)=I(j)+dt*(J(j+1)-AA(j));
%     R(j+1)=R(j)+dt*(p1*A1*inte1(j)+p1*AA1*inte11(j)-b*R(j));
    R(j+1)=R(j)+dt*(1-b)*(p1*A1*inte1(j)+p1*AA1*inte11(j));
    D(j+1)=D(j)+dt*(p2*A2*inte2(j)+p2*AA2*inte22(j));
    CJ(j+1)=CJ(j)+J(j);
  
   
end

end



for m=gap:gap
for j=1+T1*(m-1):TTT
      t(j+1)=t(j)+dt;
    CC=CC+1
  %  t(j+1)=t(j)+dt;
    h=dt;
    
    %%%%%%%%%%%beta integration%%%%%%%%%%

binte1(j)=(dt/2)*(beta(1)*J(j)+ beta(j)*J(1));
    for k=2 :j-1
        Y1(k)=dt*beta(k)*J(j-k+1);
    binte1(j)=binte1(j) + Y1(k);
    end
    
    BETA(j)=binte1(j);
    
    
    %%%%%%%%%%%recovery integration%%%%%%%%%%
    inte1(j)=(dt/2)*((exp(-t(j)/b1))*((t(j))^(a1-1))*J(1)+0);
    for k=2 :j-1
        Y1(k)=dt*(J(k))*(exp(-t(j-k+1)/b1))*(t(j-k+1))^(a1-1);
    inte1(j)=inte1(j) + Y1(k);
    end
    
   inte11(j)=(dt/2)*((exp(-t(j)/bb1))*((t(j))^(aa1-1))*J(1)+0);
    for k=2 :j-1
        Y11(k)=dt*(J(k))*(exp(-t(j-k+1)/bb1))*(t(j-k+1))^(aa1-1);
    inte11(j)=inte11(j) + Y11(k);
    end
    
    %%%%%%%%%%%dead integration%%%%%%%%%%
    inte2(j)=(dt/2)*((exp(-t(j)/b2))*(t(j))^(a2-1)*J(1)+0);
    for k=2 :j-1
      Y2(k)=dt*(J(k))*(exp(-t(j-k+1)/b2))*(t(j-k+1))^(a2-1);
    inte2(j)=inte2(j) + Y2(k);
    end
    
    inte22(j)=(dt/2)*((exp(-t(j)/bb2))*(t(j))^(aa2-1)*J(1)+0);
    for k=2 :j-1
      Y22(k)=dt*(J(k))*(exp(-t(j-k+1)/bb2))*(t(j-k+1))^(aa2-1);
    inte22(j)=inte22(j) + Y22(k);
    end
    
    
    AA(j)=p1*A1*inte1(j)+p1*AA1*inte11(j)+p2*A2*inte2(j)+p2*AA2*inte22(j);
    
    %%%%%%%%%%%level of immunity integration%%%%%%%%%%


    
    
      Ms(j)=(dt/2)*(psi(j)*(p1*A1*inte1(1)+p1*AA1*inte11(1))+psi(1)*(p1*A1*inte1(j)+p1*AA1*inte11(j)));%+2*(sum(psi(i-1:-1:2).*recovery_data(2:1:i-1))));
          
           for k=2:j-1
    
                 Ms(j)=Ms(j)+dt*psi(j-k+1)*(p1*A1*inte1(k)+p1*AA1*inte11(k));
           end
            M1(j)=(1-b)*(Ms(j)/N);
           M2(j)=eps*b*(Ms(j)/N);
           MV1(j)=0;
           MV2(j)=0;
           for n=1:m
           M1(j)=M1(j)+MM1(j-T1*(n-1));
           M2(j)=M2(j)+MM2(j-T1*(n-1));
            MV1(j)=MV1(j)+MM1(j-T1*(n-1));
           MV2(j)=MV2(j)+MM2(j-T1*(n-1));
           end
           M(j)=M1(j)+M2(j);
            MV(j)=MV1(j)+MV2(j);
           RR(j)=(dt/2)*(R(1)+R(j)+2*sum(R(2:j-1)));
           AAA(j)=dt*(p1*A1*inte1(j)+p1*AA1*inte11(j));
             RRR(j)=(dt/2)*(AAA(1)+AAA(j)+2*sum(AAA(2:j-1)));
           JJ2(j)=(dt/2)*(J2(1)+J2(j)+2*sum(J2(2:j-1)));
   % Mf(j)=MM1(j)+(Ms(j)/N);  
%      S(j+1)=S(j)+dt*(-BETA(j)*(1-Mf(j)-((I(j)+D(j))/N)));
    S(j+1)=N-(I(j)+D(j)+P(j)+M1(j)*N);
    P(j+1)=b*RRR(j)-JJ2(j)-M2(j)*N;
    J1(j+1)=(BETA(j)*S(j)/N);
    J2(j+1)=(kappa*BETA(j)*P(j)/N);
    J(j+1)=J1(j+1)+J2(j+1);
    I(j+1)=I(j)+dt*(J(j+1)-AA(j));
%     R(j+1)=R(j)+dt*(p1*A1*inte1(j)+p1*AA1*inte11(j)-b*R(j));
    R(j+1)=R(j)+dt*(1-b)*(p1*A1*inte1(j)+p1*AA1*inte11(j));
    D(j+1)=D(j)+dt*(p2*A2*inte2(j)+p2*AA2*inte22(j));
    CJ(j+1)=CJ(j)+J(j);
  
   
end

end

ICOST(:,kk)=I;
end
% 
% plot(t,I,'r','LineWidth',2);hold on

 %plot(t(1:end),I,'r','LineWidth',2);hold on


 %%%%%%%%%%% COSY FUNCTION %%%%%%%%%%%%%%
c1=0.01;c2=5;
for kk=1:length(ngap)
%  COST(kk)= (c1*sum(ICOST(:,kk))+ c2*ngap(kk));
 COST(kk)= (c1*sum(ICOST(:,kk))+ c2*ngap1(kk));
end

plot(ngap,COST,'r','LineWidth',2);hold on
plot(ngap,COST,'b.','MarkerSize',15);hold on
xlabel('Gap between two consecutive vaccinations $a$ (in months)' ,'interpreter','latex');
ylabel('Cost function $\mathcal{J}(a)$','interpreter','latex');
axis([3 30 7.3*10^5 1.6*10^6]);
