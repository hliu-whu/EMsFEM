function Dis = Linear_FEM
% reference solution:
% linear elastic analysis by using the traditional FEM on the fine-scale mesh
tic;

    Elems = load('../data/Fine_Elements.dat');
    Nodes = load('../data/Fine_Nodes.dat');
    Mod   = load('../data/Fine_Modulus.dat');
    DisBC = load('../data/Fine_DisplacementBC.dat');
    ForBC = load('../data/Fine_ForceBC.dat');
    
    emu0 = 0.3;
    thick = 1.0;
    
    bc_n = size(DisBC,1);
    bcdofs = zeros(1,bc_n);
    ld_n = size(ForBC,1);
    lddofs = zeros(ld_n,2);
    for i = 1:bc_n
        bcdofs(1,i) = (DisBC(i,1)-1)*2 + DisBC(i,2);
    end
    for i = 1:ld_n
        lddofs(i,1) = (ForBC(i,1)-1)*2 + ForBC(i,2);
        lddofs(i,2) = ForBC(i,3);
    end
    nn = size(Nodes,1);
    alldofs = 1:2*nn;
    freedofs = setdiff(alldofs,bcdofs);
    
%     lddofs(:,2) = -lddofs(:,2); 

    F = sparse(2*nn,1);
    F(lddofs(:,1),1) = lddofs(:,2);

    Dis = zeros(2*nn,1);
    

    NE = size(Elems,1);

    Kes_Elems = zeros(64,NE);
    
    for ie = 1:NE
        nods = Elems(ie,:);
        XX = Nodes(nods,1);
        YY = Nodes(nods,2);

        D = Get_D(Mod(ie,1),emu0);
        XY_s = [XX YY];
        Kes = StiffnessMatrix_FineElement(XY_s,D,thick);
        
        Kes_Elems(:,ie) = Kes(:);
    end
    
    edofMat = kron(Elems,[2,2])+repmat([-1,0],NE,4);
    iK = reshape(kron(edofMat,ones(8,1))',64*NE,1);
    jK = reshape(kron(edofMat,ones(1,8))',64*NE,1);   
    sK = reshape(Kes_Elems,64*NE,1);
    Kt = sparse(iK,jK,sK); 
    Kt = (Kt+Kt')/2;
    
    Dis(freedofs,1) = Kt(freedofs,freedofs) \ F(freedofs,1);

    OutPut(Nodes,Elems,Dis(1:2:end,1),Dis(2:2:end,1),'FEM_Linear.dat');
    
toc;
end