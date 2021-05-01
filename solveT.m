function T = solveT(u,T_start,IEN,c1,c2,kappa,p,n,nPoints,h,Tw,qw,NewtonMaxIter,Newton_tol)
%% solve Temp. field use Newton's iteration based on velocity field
T = T_start;
T(nPoints) = Tw; %% satisfy boundary
nElements = nPoints-1;
totalDOF = nPoints;
[~,rhs_s] = globalKF(IEN,c1,c2,kappa,p,n,T,u,totalDOF,nElements,h,qw);
for iter = 1:NewtonMaxIter
        
    [lhs,rhs] = globalKF(IEN,c1,c2,kappa,p,n,T,u,totalDOF,nElements,h,qw);
    dT = lhs\rhs;
    T = T+dT;
    
    %% check convergence
    [ev_rhs,ev] = Newton_check(dT,rhs,T_start,rhs_s);
    %rhs
    convergence = [ev_rhs,ev];
    
    if ev_rhs< Newton_tol(1) && ev<Newton_tol(2) 
        break
    end
    
end
end

function [ev_rhs,ev] = Newton_check(du,rhs,u_start,rhs_start)
ev_rhs = norm(rhs)/norm(rhs_start);
ev = norm(du)/norm(u_start);
end

function [lhs,rhs] = globalKF(IEN,c1,c2,kappa,p,n,T,u,totalDOF,nElements,h,qw)
nPoints = nElements+1;
lhs = zeros(totalDOF,totalDOF);
rhs = zeros(totalDOF,1);

for element = 1:nElements
    J = h(element);
    InvJ = 1/J;
    Tele = zeros(1,p+1);
    uele = zeros(1,p+1);
    for i = 1:p+1
        Tele(i) = T(IEN(element,i));
        uele(i) = u(IEN(element,i));
    end
    [k,f] = ElementStiff(c1,c2,kappa,p,n,Tele,uele,J,InvJ);
    
    %assemble    
    for i = 1:p+1
        for j = 1:p+1
            lhs(IEN(element,i),IEN(element,j)) = k(i,j)+lhs(IEN(element,i),IEN(element,j));
        end
        rhs(IEN(element,i)) = f(i)+rhs(IEN(element,i));
    end
end

%% add flux bc at y=0
rhs(1) = rhs(1)+qw;
%% apply Strong boundary condition at y=H
rhs(nPoints) = 0;
lhs(nPoints,:) = 0;
lhs(:,nPoints) = 0;
lhs(nPoints,nPoints) = 1;

end

function [k,rhs] = ElementStiff(c1,c2,kappa,p,n,T,u,J,InvJ)
nInt = 2;
[xi,w] = GaussQuad(nInt,1);
[ShapeFunc,DivSF] = ShapeFuncDiv(p,xi);
k = zeros(p+1,p+1);
f = zeros(p+1,1);

if n>1
    for i = 1:length(w)
        Nax = DivSF(:,i)*InvJ;
        Nbx = Nax';
        Na = ShapeFunc(:,i);
        ux = Nax'*u';
        Tquad = Na'*T';
        Tx = Nax'*T';
        uquad = Na'*u';
        
        k = k+Na*c2*ux^2*n*Tquad^(n-1)*J*w(i)+Nax*kappa*Nbx*J*w(i);
        f = f+Na*c2*ux^2*Tquad^n*J*w(i)+Nax*kappa*Tx*J*w(i)+Na*c1*uquad*J*w(i);
        
    end
else %% if n=0, T^n=1
    for i = 1:length(w)
        Nax = DivSF(:,i)*InvJ;
        Nbx = Nax';
        Na = ShapeFunc(:,i);
        ux = Nax'*u';
        Tquad = Na'*T';
        Tx = Nax'*T';
        uquad = Na'*u';
        
        k = k+Nax*kappa*Nbx*J*w(i);
        f = f+Na*c2*ux^2*J*w(i)+Nax*kappa*Tx*J*w(i)+Na*c1*uquad*J*w(i);
        
    end
end
rhs = -f;
end