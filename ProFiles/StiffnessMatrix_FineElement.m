function Kes = StiffnessMatrix_FineElement(XY_s,D,thick)
% Calculate the element stiffness matrix

    c=1/sqrt(3);   % integration points
    XNI=[-c,-c;c,-c;c,c;-c,c];

    Kes = zeros(8,8);       % element stiffness matrix

    for IP=1:4           % For all Gaussian points
        % form the Jacobian matrix J
        [TJ]    = Jacobian(XNI(IP,1),XNI(IP,2),XY_s);
        ITJ     = inv(TJ);
        DJ      = det(TJ); % Matrix determinant
        [pnxy]  = PNXY(ITJ,XNI(IP,1),XNI(IP,2));
        [B]     = getB(pnxy);
        ccc     = B'*D*B*DJ*thick;
        Kes     = ccc + Kes;
    end

end