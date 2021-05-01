function u= solveU(T,c1,c2,n,nPoints,h,BCval,IEN,p)
%% linearly solve velocity field based on Temp. distribution
nElement = nPoints-1;
totalDOF = nPoints;
K = zeros(totalDOF,totalDOF);
F = zeros(totalDOF,1);
for element = 1:nElement
    J = h(element);
    InvJ = 1/J;
    Tnode = [T(IEN(element,1)),T(IEN(element,end))];
    [k,f] = ElementStiff(J,InvJ,c1,c2,Tnode,n,p);
    %assemble    
    for i = 1:p+1
        for j = 1:p+1
            K(IEN(element,i),IEN(element,j)) = k(i,j)+K(IEN(element,i),IEN(element,j));
        end
        F(IEN(element,i)) = f(i)+F(IEN(element,i));
    end
    
end

%% apply Strong boundary condition
for i = 1:totalDOF
    F(i) = F(i)-K(i,1)*BCval(1)-K(i,nPoints)*BCval(2);
end
F(1) = BCval(1);
F(nPoints) = BCval(2);

K(1,:) = 0;
K(:,1) = 0;
K(1,1) = 1;
K(nPoints,:) = 0;
K(:,nPoints) = 0;
K(nPoints,nPoints) = 1;
%conditionK = cond(K);
%K

u = K\F;
end

function [k,f] = ElementStiff(J,InvJ,c1,c2,Tnode,n,p)
%p=1;%% use linear element
nInt = 2; %% use 2 intergration points
[xi,w] = GaussQuad(nInt,1);
[ShapeFunc,DivSF] = ShapeFuncDiv(p,xi);
k = zeros(p+1,p+1);
f = zeros(p+1,1);

for i = 1:length(w)
    Nax = DivSF(:,i)*InvJ;
    Nbx = Nax';
    Na = ShapeFunc(:,i);
    T = Na'*Tnode';
    coef = c2*T^n;

    k = k +Nax*coef*Nbx*J*w(i);
    f = f + Na*c1*J*w(i);
end
end