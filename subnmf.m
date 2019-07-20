function [H,grad,iter,r]=subnmf(V,W,H0,tol,iterlimit,deltal0,deltal1,rho,yita,tr1)
ep1=10^(-7);
WtV = W'*V; WtW = W'*W;
H = H0;
gradh0 = WtW*H - WtV;%Gradient matrix
tr2=sum(sum(H.*WtV));
tr3=sum(sum(WtW.*(H*H')));
f0=(tr1-2*tr2+tr3)/2;
d0=-gradh0;%Direction matrix
pd=norm(gradh0);
for t=1:iterlimit
    % stopping condition
     if pd < tol 
         break
     end
    [H1,f1]=NMF_armijol(deltal0,deltal1,rho,WtW,WtV,H0,d0,f0,gradh0,tr1);
    x0=sum(sum(H0.^2));
    x1=sum(sum(H1.^2));
    aa=min(x0,x1);
    bb=min(f0,f1);
    %% avoid infinite loop
    if abs(f0-f1)/bb<ep1 && abs(x0-x1)/aa<ep1
        r=1;
        break;
    end 
    gradh1 = WtW*H1 - WtV;
    yk_1=gradh1-gradh0;
    sk_1=H1-H0; %matrix
    dk_1=d0;
    dyk_1=sum(sum(dk_1.*yk_1));
    d_1g1=sum(sum(dk_1.*gradh1));
    g1y_1=sum(sum((gradh1.*yk_1)));
    g1s_1=sum(sum(gradh1.*sk_1));
    ngradh1=sum(sum(gradh1.^2));
    ngradh0=sum(sum(gradh0.^2));
    dk_1yk_1=dyk_1-d_1g1*g1y_1/ngradh1;
    g1g0=sum(sum(gradh1.*gradh0));
    if dk_1yk_1>yita*ngradh0
        sitak=(dyk_1-(d_1g1)*g1s_1/ngradh1)/dk_1yk_1;
        betak=(g1y_1-g1s_1)/dk_1yk_1;
    else
        sitak=dyk_1/ngradh0-(d_1g1)*(g1g0)/ngradh1/ngradh0;
        betak=g1y_1/ngradh0;
    end
    d1=-sitak*gradh1+betak*dk_1;
    H0=H1;
    d0=d1;
    f0=f1;
    gradh0 =gradh1;
    pd=norm(gradh0);
end
iter=t;
grad=gradh0;
H=H0;
r=0;
end