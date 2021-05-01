function [ShapeFunc,DivSF,SF2,SFP] = ShapeFuncDiv(p,x)

[ShapeFunc,DivSF,SF2] = ShapeFuncDivBernstein(p,x);
if p==2
    SFP=[2;-4;2];
elseif p==3
    SFP = [-6;18;-18;6];
end



end

function [ShapeFunc,DivSF,SF2] = ShapeFuncDivBernstein(p,x)
L = length(x);
DivSF = zeros(p+1,L);
SF2 = zeros(p+1,L);
if p==1
    DivSF(1,:) = -1;
    DivSF(2,:) = 1;
    SF2(1,:) = 0;
    SFP(2,:) = 0;
else
    for i = 1:length(x)
        xi = x(i);
        Bpn1 = [0,bernsteinMatrix(p-1,xi),0];
        Bpn2 = [0,0,bernsteinMatrix(p-2,xi),0,0];
        dB = zeros(p+1,1);
        dB2 = zeros(p+1,1);
        for k = 1:p+1
            dB(k) = p*(Bpn1(k) - Bpn1(k+1));
            dB2(k) = p*(p-1)*(Bpn2(k)-2*Bpn2(k+1)+Bpn2(k+2));
        end
        DivSF(:,i) = dB;
        SF2(:,i) = dB2;
    end
end
B = bernsteinMatrix(p,x);
ShapeFunc = B';
end