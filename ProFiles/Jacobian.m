function [TJ]=Jacobian(ex,ey,XX)
% compute the Jacobian matrix of 4-node Q4 element

    pn=zeros(2,4);
    pn(1,1)=-0.25*(1-ey);     pn(1,2)=0.25*(1-ey);
    pn(1,3)=0.25*(1+ey);      pn(1,4)=-0.25*(1+ey);
    pn(2,1)=-0.25*(1-ex);     pn(2,2)=-0.25*(1+ex);
    pn(2,3)=0.25*(1+ex);      pn(2,4)=0.25*(1-ex);

    TJ=pn*XX;

end