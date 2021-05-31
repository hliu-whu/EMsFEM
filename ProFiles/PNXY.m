function pnxy = PNXY(ITJ,ex,ey)
% Calculate the partial derivative of the shape function to the global coordinate

    pn=zeros(2,4);
    pn(1,1)=-0.25*(1-ey);    pn(1,2)=0.25*(1-ey);
    pn(1,3)=0.25*(1+ey);     pn(1,4)=-0.25*(1+ey);
    pn(2,1)=-0.25*(1-ex);    pn(2,2)=-0.25*(1+ex);
    pn(2,3)=0.25*(1+ex);     pn(2,4)=0.25*(1-ex);
    
    pnxy = ITJ*pn;
    
end