function [Dis_Coarse,Dis_Fine] = Linear_EMs_Heter(type)
% Linear elastic analysis of non-periodic structures by using the EMsFEM
tic;

    C_Elems = load('../data/Coarse_Elements.dat');
    C_Nodes = load('../data/Coarse_Nodes.dat');
    S_Elems = load('../data/Sub_Elements.dat');
    S_Nodes = load('../data/Sub_Nodes.dat');
    DisBC = load('../data/Coarse_DisplacementBC.dat');
    ForID = load('../data/Coarse_ForceID.dat');
    
    [MatMappingOS,Sub2OSElems] = MatMappingOS_fun(C_Nodes,C_Elems,1e-6);
    
    emu0 = 0.3;
    thick = 1.0;
    bc_n = size(DisBC,1);
    bcdofs = zeros(1,bc_n);
    for i = 1:bc_n
        bcdofs(1,i) = (DisBC(i,1)-1)*2 + DisBC(i,2);
    end
    
    nn = size(C_Nodes,1);
    alldofs = 1:2*nn;
    freedofs = setdiff(alldofs,bcdofs);
    
    Dis_Coarse = zeros(2*nn,1);
    
    NE = size(C_Elems,1);
    dofs = zeros(8,1);
    SN_Cells = cell(NE,1);
    CEKs_Cells = zeros(64,NE);
    CEFs_Cells = zeros(8,NE);

    for ie = 1:NE
        nods = C_Elems(ie,:);
        xdofs = 2*nods-1;
        ydofs = 2*nods;
        dofs(1:2:8,1) = xdofs;
        dofs(2:2:8,1) = ydofs;
        OE = MatMappingOS(:,ie);
        SE = OE(Sub2OSElems,1);
        if(type==1)
            [SN,CEKs,CEFs] = GetShapeFunsOS4(S_Nodes,S_Elems,SE,OE,1e-5,1e-5,emu0,thick);
        else
            [SN,CEKs,CEFs] = GetShapeFunsOP4(S_Nodes,S_Elems,SE,OE,1e-5,1e-5,emu0,thick);
        end
        SN_Cells{ie,1} = SN;
        CEKs_Cells(:,ie) = CEKs(:);
        CEFs_Cells(:,ie) = CEFs(:);
    end

    edofMat_cell = kron(C_Elems,[2,2])+repmat([-1,0],NE,4);
    iK_cell = reshape(kron(edofMat_cell,ones(8,1))',64*NE,1);
    jK_cell = reshape(kron(edofMat_cell,ones(1,8))',64*NE,1);   
    sK_cell = reshape(CEKs_Cells,64*NE,1);
    Kt = sparse(iK_cell,jK_cell,sK_cell); 
    Kt = (Kt+Kt')/2;
    
    iF_cell = reshape(edofMat_cell',8*NE,1);
    sF_cell = reshape(CEFs_Cells.*(ones(8,1)*ForID'),8*NE,1);
    F = sparse(iF_cell,1,sF_cell); 
%     F = F*8e3; 

    
    Dis_Coarse(freedofs,1) = Kt(freedofs,freedofs) \ F(freedofs,1);


    OutPut(C_Nodes,C_Elems,Dis_Coarse(1:2:end,1),Dis_Coarse(2:2:end,1),'EMs_Coarse_Linear.dat');
     
    Dis_Fine = GetFineRes_Hetero(C_Nodes,C_Elems,S_Nodes,Dis_Coarse,SN_Cells);
    
toc;
end