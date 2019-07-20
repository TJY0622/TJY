clc
clear
[AA,BB]=xlsread('E:\2019_7_8 - Copy.xlsx') ;
A=AA;
% A=[mean(AA(:,1:3),2),mean(AA(:,4:6),2),mean(AA(:,7:9),2),mean(AA(:,10:12),2),...
%       mean(AA(:,13:15),2),mean(AA(:,16:18),2),mean(AA(:,19:21),2),mean(AA(:,22:24),2)];
% aa=2;
% A=[AA(:,aa),AA(:,aa+3),AA(:,aa+6),AA(:,aa+9),AA(:,aa+12),AA(:,aa+15),AA(:,aa+18),AA(:,aa+21)];
% aa=1;
% A=[AA(:,aa),AA(:,aa+3),AA(:,aa+6),AA(:,aa+9),AA(:,aa+12),AA(:,17),AA(:,20),AA(:,23)];
% A=[AA(:,2:3),AA(:,5:6),AA(:,8:9),AA(:,11:12),...
%       AA(:,14:15),AA(:,17:18),AA(:,20:21),AA(:,23:24)];
% A=[AA(:,1:2),AA(:,4:5),AA(:,7:8),AA(:,10:11),...
%       AA(:,13:14),AA(:,16:17),AA(:,19:20),AA(:,22:23)];
% A=[AA(:,1),AA(:,3),AA(:,4),AA(:,6),AA(:,7),AA(:,9),AA(:,10),AA(:,12),...
%     AA(:,13),AA(:,15),AA(:,16),AA(:,18),AA(:,19),AA(:,21),AA(:,22),AA(:,24)];
%%exclude all rows that are 0
[n1,l1]=size(A);
%%standardization
b=max(A,[],2);
a=min(A,[],2);
V=(A-kron(a,ones(1,l1)))./kron((b-a),ones(1,l1));
%%exclude all rows that are 0
[nn,mm]=find(isnan(V)==1);
V(nn,:)=[];
BB(nn+1,:)=[];
[nnn,lll]=size(V);
s=6;
n=nnn;m=s;l=lll;
%%initialzation
W=abs(randn(nnn,s));
H=abs(randn(s,lll));
%%parameter setting
iterlimit=10^4;
maxiter=10^3;
timelimit=2*10^3;
tol=10^-11;
epsilon=10^-7;
% deltal0W=0.4;deltal1W=0.0001;rhoW=0.65;yitaW=0.001;
deltal0=0.4;
deltal1=10^-4;
rho=0.65;
yita=10^-3;
[iter,totaliter, projnorm,time,W0,H0]=nmf(V,W,H,iterlimit,maxiter,timelimit,tol,epsilon,deltal0,deltal1,rho,yita);
%% the composition matrix R
record=zeros(s,lll);
[nnn1,lll1]=size(W0);
for i=1:lll
	HH=repmat(H0(:,i)',nnn1,1);
	first=W0.*HH;
	sum1=sum(first,2);
    n=find(sum1==0);
    first(n,:)=[];
    sum1(n,:)=[];
	D1=first./repmat(sum1,1,lll1);
    nn=length(D1);
	P1=sum(D1,1)/nn;
	P1=P1';
	record(:,i)=P1;
end