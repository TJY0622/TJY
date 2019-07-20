function [H,f]=NMF_armijol(deltal0,deltal1,rho,WtW,WtV,H0,d0,f0,gradh0,tr1)
alpha=1;
% [n,m]=size(W0);
% D0=reshape(d0,n,m); 
H1=H0+alpha*d0;
% tr2=sum(sum(H.*WtV));
% tr3=sum(sum(WtW.*(H*H')));
% f0=(-2*tr2+tr3)/2;
tr2=sum(sum(H1.*WtV));
tr3=sum(sum(WtW.*(H1*H1')));
fsub1=tr1-2*tr2+tr3;
f1=fsub1/2;
gd=sum(sum(gradh0.*d0));
tr4=sum(sum(d0.^2));
N=0;
 fmpprp=f0+deltal0*alpha*gd-deltal1*(alpha^2)*tr4;
 while f1>fmpprp
     if N > 100
         break;
     end
     alpha=rho*alpha;
     H1=H0+alpha*d0;   %matrix
     fmpprp=f0+deltal0*alpha*gd-deltal1*(alpha^2)*tr4;
     trr2=sum(sum(H1.*WtV));
     trr3=sum(sum(WtW.*(H1*H1')));
     fsub1=tr1-2*trr2+trr3;
     f1=fsub1/2;
     N=N+1;
 end
 H=H1;
 f=f1;
end
