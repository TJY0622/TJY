function [iter,totaliter, projnorm,time,W,H]=nmf(V,W,H,iterlimit,maxiter,timelimit,tol,epsilon,deltal0,deltal1,rho,yita)

%%
% resid0=norm(V-W*H,'fro')/norm(V,'fro');
%W = Winit; H = Hinit; initt = cputime;
totaliter=0;
gradW = W*(H*H') - V*H'; gradH = (W'*W)*H - W'*V;
tr1=sum(sum(V.^2));
initgrad = norm([gradW; gradH'],'fro');
fprintf('Init gradient norm %f\n', initgrad);
% tolW = max(epsilon,tol)*initgrad; 
tolW = epsilon; 
tolH = tolW;
projnorm0 = norm([gradW(gradW<0 | W>0); gradH(gradH<0 | H>0)]); %projected gradient

%%
tic
initt=cputime;
for iter=1:maxiter
    % stopping condition
    if (projnorm0 < tol*initgrad) || ((cputime-initt) > timelimit)
        break;
    end
    %get w first
    [W,gradW,iterW,r]=subnmf(V',H',W',tolW,iterlimit,deltal0,deltal1,rho,yita,tr1);
    W=max(W, 0);
    W = W'; gradW = gradW';   %subproblem
    if iterW<=5 && r==0
        tolW = 0.1 * tolW;
        fprintf('iterW<=5\n')
    end
    totaliter=totaliter+iterW;
    [H,gradH,iterH,r]=subnmf(V,W,H,tolH,iterlimit,deltal0,deltal1,rho,yita,tr1);
    H=max(H, 0);
    if iterH<=5 && r==0
        tolH = 0.1 * tolH;
         fprintf('iterH<=5\n')
    end
    totaliter=totaliter+iterH;
    projnorm1 = norm([gradW(gradW<0 | W>0); gradH(gradH<0 | H>0)]); %projected gradient
    projnorm0=projnorm1;
end
time=toc;
projnorm=projnorm0;
end