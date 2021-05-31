function [Dis_Coarse,Dis_Fine] = Linear_EMs
% Linear elastic analysis of periodic structures by using the EMsFEM
tic;
    
    C_Elems = load('../data/Coarse_Elements.dat');
    C_Nodes = load('../data/Coarse_Nodes.dat');
    S_Elems = load('../data/Sub_Elements.dat');
    S_Nodes = load('../data/Sub_Nodes.dat');
    SE = load('../data/Sub_Modulus.dat');
    
    emu0 = 0.3;
    thick = 1.0;
    CNE = size(C_Elems,1);
    edofMat_coarsegrid = kron(C_Elems,[2,2])+repmat([-1,0],CNE,4);
    iK_coarsegrid = reshape(kron(edofMat_coarsegrid,ones(8,1))',64*CNE,1);
    jK_coarsegrid = reshape(kron(edofMat_coarsegrid,ones(1,8))',64*CNE,1);   
    SNE = size(S_Elems,1);
    edofMat_subgrid = kron(S_Elems,[2,2])+repmat([-1,0],SNE,4);    
    iK_subgrid = reshape(kron(edofMat_subgrid,ones(8,1))',64*SNE,1);
    jK_subgrid = reshape(kron(edofMat_subgrid,ones(1,8))',64*SNE,1); 
    
    EK_subgrid = zeros(64,SNE);
    
    for ise = 1:SNE
        E = SE(ise);
        D = Get_D(E,emu0);
        EK = StiffnessMatrix_FineElement(S_Nodes(S_Elems(ise,:),:),D,thick);
        EK_subgrid(:,ise) = EK(:);
    end
    K_subgrid = sparse(iK_subgrid,jK_subgrid,EK_subgrid(:)); 
    K_subgrid = (K_subgrid+K_subgrid')/2;
    
    [SN,CEKs,CEFs] = GetShapeFunsP4(S_Nodes,S_Elems,K_subgrid,1e-5,1e-5); % construct the numerical shape function with the periodic boundary condition
    
    sK_coarsegrid = reshape(CEKs(:)*ones(1,CNE),64*CNE,1);
    K_coarsegrid = sparse(iK_coarsegrid,jK_coarsegrid,sK_coarsegrid);
    K_coarsegrid = (K_coarsegrid+K_coarsegrid')/2;
    
    DisBC = load('../data/Coarse_DisplacementBC.dat');
    ForID = load('../data/Coarse_ForceID.dat');
    bc_n = size(DisBC,1);
    bcdofs = zeros(1,bc_n);
    for i = 1:bc_n
        bcdofs(1,i) = (DisBC(i,1)-1)*2 + DisBC(i,2);
    end
    nn = size(C_Nodes,1);
    alldofs = 1:2*nn;
    freedofs = setdiff(alldofs,bcdofs);
    CNE = size(C_Elems,1);
    iF_cell = reshape(edofMat_coarsegrid',8*CNE,1);
    sF_cell = reshape((CEFs*ForID'),8*CNE,1);
    F = sparse(iF_cell,1,sF_cell); 
    
    Dis_Coarse = zeros(2*nn,1);
    Dis_Coarse(freedofs,1) = K_coarsegrid(freedofs,freedofs) \ F(freedofs,1);

    OutPut(C_Nodes,C_Elems,Dis_Coarse(1:2:end,1),Dis_Coarse(2:2:end,1),'EMs_Coarse_Linear.dat');
    Dis_Fine = GetFineRes(C_Nodes,C_Elems,S_Nodes,Dis_Coarse,SN);
    
toc;
end