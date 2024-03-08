

% bar(days,AA,'r','LineWidth',1);
%  hold on
%  

 
function vaccine_induced_immunity
tic
clear all
 close all
% 
% format long


set(0,'DefaultAxesFontSize',11);

%load('vaccine_immunity_1.mat');
BB=[10 40 61.69 67.9 76 122.6 128.7 135.8 176.35 188.5 215.9 213.85 226 217.9 225 234.12 241.26 240.2 224 226 242 244.4 258.4 191.5 138.9 115.5 ...
           110 110 100 105 105 100 90 80 80 50 30 30 30 10 10 5 ...
];
A=BB/max(BB);
days= 1:7:294;
%AA=average/max(average);
x=1:7:294;%0:0.01:5;
%a=2;b=3;
B=0*x;
b=2;d=3;
%A=((b*d)/(d-b))*(exp(-b*x)-exp(-d*x)); %(a*x)./(1+b*x);
%                                      
% bar(x,A,'b','LineWidth',1);hold on
plot(x,A,'b.','LineWidth',25);hold on
%plot(x,A,'b.','MarkerSize',25);hold on
a1=100;b1=1;c1=10;%INITIAL GUESS OF a1, b1
function error_in_data = moder(k)
    a1 = k(1);
    b1 = k(2);
    c1 = k(3);
    x=1:7:294;
    for l=1:length(x)
    B(l) = a1*exp((-(x(l)-b1).^2)./(2*c1*c1));%c1*((b1*a1)/(a1-b1))*(exp(-b1*x(l))-exp(-a1*x(l)));
    end
error_in_data = sum((A - B).^2); %computes SSE
end
 k = [a1 b1 c1]; 
[k,fval] = fminsearch(@moder,k); 
a1
b1
c1
xx=1:1:294;
% plot(x,B,'r','LineWidth',2);hold on
plot(xx,a1*exp((-(xx-b1).^2)./(2*c1*c1)),'r','LineWidth',3);hold on
axis([0 294 0 1]);
xlabel('Days post vaccination $\rightarrow$','interpreter','latex');
ylabel('effectiveness of vaccine-induced immunity $\rightarrow$','interpreter','latex');
% xt={'15/11/2021'  ; '17/1/2022' ; '21/3/2022'; '23/5/2022'} ; 
% set(gca,'xtick',[14 85 148 213]); 
% set(gca,'xticklabel',xt);
toc
end