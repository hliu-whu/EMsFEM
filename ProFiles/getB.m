function B0 = getB(pnxy)
% calculate the strain-displacement matrix

    B0=zeros(3,8);
    for i=1:4
        B0(1,2*i-1) = pnxy(1,i);
        B0(1,2*i)   = 0;
        B0(2,2*i-1) = 0;
        B0(2,2*i)   = pnxy(2,i);
        B0(3,2*i-1) = pnxy(2,i);
        B0(3,2*i)   = pnxy(1,i);
    end
    
end