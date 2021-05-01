function run()
c1 = 1;
c2 = 0.1;
kappa = 1;
n = 2;
H=1;
dx = 0.1;
Tw = 1;
qw = 1;
BCvalu = [0,0];
nIter = 5;
NewtonMaxIter = 10;
Newton_tol = [1e-4,1e-4];
[u,T]=solveChannel(c1,c2,kappa,n,dx,H,Tw,qw,BCvalu,nIter,NewtonMaxIter,Newton_tol);
figure(1)
y = 0:dx:H;
hold on
plot(u,y,'b-','LineWidth',2)
hold on
plot(T,y,'r-','LineWidth',2)

%Tverify = verifyT(y);
%hold on 
%plot(Tverify,y,'k')
end

function [u,T]=solveChannel(c1,c2,kappa,n,dx,H,Tw,qw,BCvalu,nIter,NewtonMaxIter,Newton_tol)
nElements = H/dx;
nGrid = nElements+1;
T = ones(nGrid,1)*Tw;
h = ones(nGrid-1,1)*dx;
p=1;
IEN = getIEN(nElements,nGrid,p); %% connectivity array
%% solve U and T iteratively
for iter = 1:nIter
    u = solveU(T,c1,c2,n,nGrid,h,BCvalu,IEN,p); %% update u based on current T
    Tnew = solveT(u,T,IEN,c1,c2,kappa,p,n,nGrid,h,Tw,qw,NewtonMaxIter,Newton_tol); %% update T based on current U 
    unew = solveU(Tnew,c1,c2,n,nGrid,h,BCvalu,IEN,p);
    [u_conv(iter),T_conv(iter)] = convergence_check(u,unew,T,Tnew); %% check convergence
    T=Tnew;
end
figure(2)
loglog (1:nIter,u_conv,'b-','LineWidth',2)
hold on
loglog(1:nIter,T_conv,'r-','LineWidth',2)
end


function IEN = getIEN(nElements,nPoints,p)
IEN = zeros(nElements,p+1);
for ele = 1:nElements
    IEN(ele,1) = ele;
    IEN(ele,end) = ele+1;
    IEN(ele,2:p) = nPoints+(ele-1)*(p-1)+1:nPoints+(ele-1)*(p-1)+p-1;
end
end

function [u_conv,T_conv] = convergence_check(u,unew,T,Tnew)
u_conv = norm(u-unew)/norm(u);
T_conv = norm(T-Tnew)/norm(T);
end

function T = verifyT(y)
T = (1/12)*y.^3-(1/8)*y.^2+y+0.0825;
end
