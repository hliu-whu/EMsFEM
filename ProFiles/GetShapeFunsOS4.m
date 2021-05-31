function [SN,CEKs,CEFs] = GetShapeFunsOS4(S_Nodes,S_Elems,SE,OE,sx,sy,mu,t)

% Load the node information of the oversampling mesh and the sub-grid mesh
%==========================================================================
    O_Nodes = load('..\data\OS_Nodes.dat');
    O_Elems = load('..\data\OS_Elements.dat');
    
    C = 1.0e6;   
    
    Xmin = min(O_Nodes(:,1)); % The coordinates of the four corners
    Ymin = min(O_Nodes(:,2));
    Xmax = max(O_Nodes(:,1));
    Ymax = max(O_Nodes(:,2));
    edge1 = find(abs(O_Nodes(:,2)-Ymin)<sy);
    edge2 = find(abs(O_Nodes(:,1)-Xmax)<sx);
    edge3 = find(abs(O_Nodes(:,2)-Ymax)<sy);
    edge4 = find(abs(O_Nodes(:,1)-Xmin)<sx);
    side = [edge1(:,1) O_Nodes(edge1(:,1),1)];
    side = sortrows(side,2);
    s1 = side(:,1);              
    side = [edge2(:,1) O_Nodes(edge2(:,1),2)];
    side = sortrows(side,2);
    s2 = side(:,1);              
    side = [edge3(:,1) O_Nodes(edge3(:,1),1)];
    side = sortrows(side,-2);
    s3 = side(:,1);              
    side = [edge4(:,1) O_Nodes(edge4(:,1),2)];
    side = sortrows(side,-2);
    s4 = side(:,1);              
    
    ONN = size(O_Nodes,1);
    ONE = size(O_Elems,1);
    SNN = size(S_Nodes,1);
    SNE = size(S_Elems,1);
    OF = sparse(2*ONN,1);
        
    EK_OS_Elems = zeros(64,ONE);
    for ioe = 1:ONE
        E = OE(ioe);
        D = Get_D(E,mu);
        EK = StiffnessMatrix_FineElement(O_Nodes(O_Elems(ioe,:),:),D,t);
        EK_OS_Elems(:,ioe) = EK(:);
    end
    edofMat_os = kron(O_Elems,[2,2])+repmat([-1,0],ONE,4);
    iK_os = reshape(kron(edofMat_os,ones(8,1))',64*ONE,1);
    jK_os = reshape(kron(edofMat_os,ones(1,8))',64*ONE,1);
    sK_os = reshape(EK_OS_Elems,64*ONE,1);
    OK = sparse(iK_os,jK_os,sK_os);
    OK = (OK+OK')/2;
    
    TOK = OK;
    TOF = OF;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % construction of the temporary shape function of the oversampling mesh
    % u1
    for i = 2:size(s2,1)
        TOK(s2(i,1)*2-1,:) = 0;
        TOK(:,s2(i,1)*2-1) = 0;
        TOK(s2(i,1)*2-1,s2(i,1)*2-1) = 1;
        TOF(s2(i,1)*2-1,1) = 0;
    end
    for i = 2:size(s3,1)
        TOK(s3(i,1)*2-1,:) = 0;
        TOK(:,s3(i,1)*2-1) = 0;
        TOK(s3(i,1)*2-1,s3(i,1)*2-1) = 1;
        TOF(s3(i,1)*2-1,1) = 0;
    end
    TOK(s3(1,1)*2,:) = 0;
    TOK(:,s3(1,1)*2) = 0;
    TOK(s3(1,1)*2,s3(1,1)*2) = 1;
    TOF(s3(1,1)*2,1) = 0;
    for i = 2:size(s1,1)
        TOK(s1(i,1)*2-1,s1(i,1)*2-1) = C*TOK(s1(i,1)*2-1,s1(i,1)*2-1);
        TOF(s1(i,1)*2-1,1) = TOK(s1(i,1)*2-1,s1(i,1)*2-1)*(1-(i-1)/(size(s1,1)-1));
    end
    for i = 2:size(s4,1)
        TOK(s4(i,1)*2-1,s4(i,1)*2-1) = C*TOK(s4(i,1)*2-1,s4(i,1)*2-1);
        TOF(s4(i,1)*2-1,1) = TOK(s4(i,1)*2-1,s4(i,1)*2-1)*((i-1)/(size(s4,1)-1));
    end
    U = TOK \ TOF;
    ON1x = U(1:2:end);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % u2
    TOK = OK;
    TOF = OF;
    for i = 2:size(s4,1)
        TOK(s4(i,1)*2-1,:) = 0;
        TOK(:,s4(i,1)*2-1) = 0;
        TOK(s4(i,1)*2-1,s4(i,1)*2-1) = 1;
        TOF(s4(i,1)*2-1,1) = 0;
    end
    TOK(s4(1,1)*2,:) = 0;
    TOK(:,s4(1,1)*2) = 0;
    TOK(s4(1,1)*2,s4(1,1)*2) = 1;
    TOF(s4(1,1)*2,1) = 0;
    for i = 2:size(s3,1)
        TOK(s3(i,1)*2-1,:) = 0;
        TOK(:,s3(i,1)*2-1) = 0;
        TOK(s3(i,1)*2-1,s3(i,1)*2-1) = 1;
        TOF(s3(i,1)*2-1,1) = 0;
    end
    for i = 2:size(s1,1)
        TOK(s1(i,1)*2-1,s1(i,1)*2-1) = C*TOK(s1(i,1)*2-1,s1(i,1)*2-1);
        TOF(s1(i,1)*2-1,1) = TOK(s1(i,1)*2-1,s1(i,1)*2-1)*((i-1)/(size(s1,1)-1));
    end
    for i = 2:size(s2,1)
        TOK(s2(i,1)*2-1,s2(i,1)*2-1) = C*TOK(s2(i,1)*2-1,s2(i,1)*2-1);
        TOF(s2(i,1)*2-1,1) = TOK(s2(i,1)*2-1,s2(i,1)*2-1)*(1-(i-1)/(size(s2,1)-1));
    end
    U = TOK \ TOF;
    ON2x = U(1:2:end);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % u3
    TOK = OK;
    TOF = OF;
    for i = 2:size(s4,1)
        TOK(s4(i,1)*2-1,:) = 0;
        TOK(:,s4(i,1)*2-1) = 0;
        TOK(s4(i,1)*2-1,s4(i,1)*2-1) = 1;
        TOF(s4(i,1)*2-1,1) = 0;
    end
    TOK(s1(1,1)*2,:) = 0;
    TOK(:,s1(1,1)*2) = 0;
    TOK(s1(1,1)*2,s1(1,1)*2) = 1;
    TOF(s1(1,1)*2,1) = 0;
    for i = 2:size(s1,1)
        TOK(s1(i,1)*2-1,:) = 0;
        TOK(:,s1(i,1)*2-1) = 0;
        TOK(s1(i,1)*2-1,s1(i,1)*2-1) = 1;
        TOF(s1(i,1)*2-1,1) = 0;
    end
    for i = 2:size(s3,1)
        TOK(s3(i,1)*2-1,s3(i,1)*2-1) = C*TOK(s3(i,1)*2-1,s3(i,1)*2-1);
        TOF(s3(i,1)*2-1,1) = TOK(s3(i,1)*2-1,s3(i,1)*2-1)*(1-(i-1)/(size(s3,1)-1));
    end
    for i = 2:size(s2,1)
        TOK(s2(i,1)*2-1,s2(i,1)*2-1) = C*TOK(s2(i,1)*2-1,s2(i,1)*2-1);
        TOF(s2(i,1)*2-1,1) = TOK(s2(i,1)*2-1,s2(i,1)*2-1)*((i-1)/(size(s2,1)-1));
    end
    U = TOK \ TOF;
    ON3x = U(1:2:end);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % u4
    TOK = OK;
    TOF = OF;
    for i = 2:size(s1,1)
        TOK(s1(i,1)*2-1,:) = 0;
        TOK(:,s1(i,1)*2-1) = 0;
        TOK(s1(i,1)*2-1,s1(i,1)*2-1) = 1;
        TOF(s1(i,1)*2-1,1) = 0;
    end
    TOK(s2(1,1)*2,:) = 0;
    TOK(:,s2(1,1)*2) = 0;
    TOK(s2(1,1)*2,s2(1,1)*2) = 1;
    TOF(s2(1,1)*2,1) = 0;
    for i = 2:size(s2,1)
        TOK(s2(i,1)*2-1,:) = 0;
        TOK(:,s2(i,1)*2-1) = 0;
        TOK(s2(i,1)*2-1,s2(i,1)*2-1) = 1;
        TOF(s2(i,1)*2-1,1) = 0;
    end
    for i = 2:size(s3,1)
        TOK(s3(i,1)*2-1,s3(i,1)*2-1) = C*TOK(s3(i,1)*2-1,s3(i,1)*2-1);
        TOF(s3(i,1)*2-1,1) = TOK(s3(i,1)*2-1,s3(i,1)*2-1)*((i-1)/(size(s3,1)-1));
    end
    for i = 2:size(s4,1)
        TOK(s4(i,1)*2-1,s4(i,1)*2-1) = C*TOK(s4(i,1)*2-1,s4(i,1)*2-1);
        TOF(s4(i,1)*2-1,1) = TOK(s4(i,1)*2-1,s4(i,1)*2-1)*(1-(i-1)/(size(s4,1)-1));
    end
    U = TOK \ TOF;
    ON4x = U(1:2:end);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    TS_Nodes = S_Nodes;
    Xmin = min(TS_Nodes(:,1)); 
    Ymin = min(TS_Nodes(:,2));
    Xmax = max(TS_Nodes(:,1));
    Ymax = max(TS_Nodes(:,2));
    VNos = VerticesNodesNo(TS_Nodes,Xmin,Ymin,Xmax,Ymax,sx,sy);
    % sub2os
    sub2os = zeros(SNN,1);
    for isn = 1:SNN
        snx = S_Nodes(isn,1);
        sny = S_Nodes(isn,2);
        osnum = find((abs(snx-O_Nodes(:,1))<sx) & (abs(sny-O_Nodes(:,2))<sy));
        [a,b] = size(osnum);
        if(a*b~=1)
            disp('Error in Function Get_ShapeFunNxy_OverSampling!!!');
        else
            sub2os(isn,1) = osnum;
        end
    end
    mc = zeros(4,4);
    for i = 1:4
        i1 = sub2os(VNos(i,1),1);
        mc(1,i) = ON1x(i1);
        mc(2,i) = ON2x(i1);
        mc(3,i) = ON3x(i1);
        mc(4,i) = ON4x(i1);
    end
    mc = inv(mc);
    SN1x = zeros(SNN,1);
    SN2x = zeros(SNN,1);
    SN3x = zeros(SNN,1);
    SN4x = zeros(SNN,1);
    for i = 1:SNN
        i1 = sub2os(i,1);
        SN1x(i,1) = mc(1,1)*ON1x(i1)+mc(1,2)*ON2x(i1)+mc(1,3)*ON3x(i1)+mc(1,4)*ON4x(i1);
        SN2x(i,1) = mc(2,1)*ON1x(i1)+mc(2,2)*ON2x(i1)+mc(2,3)*ON3x(i1)+mc(2,4)*ON4x(i1);
        SN3x(i,1) = mc(3,1)*ON1x(i1)+mc(3,2)*ON2x(i1)+mc(3,3)*ON3x(i1)+mc(3,4)*ON4x(i1);
        SN4x(i,1) = mc(4,1)*ON1x(i1)+mc(4,2)*ON2x(i1)+mc(4,3)*ON3x(i1)+mc(4,4)*ON4x(i1);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % v1
    TOK = OK;
    TOF = OF;
    for i = 2:size(s2,1)
        TOK(s2(i,1)*2,:) = 0;
        TOK(:,s2(i,1)*2) = 0;
        TOK(s2(i,1)*2,s2(i,1)*2) = 1;
        TOF(s2(i,1)*2,1) = 0;
    end
    TOK(s3(1,1)*2-1,:) = 0;
    TOK(:,s3(1,1)*2-1) = 0;
    TOK(s3(1,1)*2-1,s3(1,1)*2-1) = 1;
    TOF(s3(1,1)*2-1,1) = 0;
    for i = 2:size(s3,1)
        TOK(s3(i,1)*2,:) = 0;
        TOK(:,s3(i,1)*2) = 0;
        TOK(s3(i,1)*2,s3(i,1)*2) = 1;
        TOF(s3(i,1)*2,:) = 0;
    end
    for i = 2:size(s1,1)
        TOK(s1(i,1)*2,s1(i,1)*2) = C*TOK(s1(i,1)*2,s1(i,1)*2);
        TOF(s1(i,1)*2,1) = TOK(s1(i,1)*2,s1(i,1)*2)*(1-(i-1)/(size(s1,1)-1));
    end
    for i = 2:size(s4,1)
        TOK(s4(i,1)*2,s4(i,1)*2) = C*TOK(s4(i,1)*2,s4(i,1)*2);
        TOF(s4(i,1)*2,1) = TOK(s4(i,1)*2,s4(i,1)*2)*((i-1)/(size(s4,1)-1));
    end
    V = TOK \ TOF;
    ON1y = V(2:2:end);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % v2
    TOK = OK;
    TOF = OF;
    for i = 2:size(s4,1)
        TOK(s4(i,1)*2,:) = 0;
        TOK(:,s4(i,1)*2) = 0;
        TOK(s4(i,1)*2,s4(i,1)*2) = 1;
        TOF(s4(i,1)*2,1) = 0;
    end
    TOK(s4(1,1)*2-1,:) = 0;
    TOK(:,s4(1,1)*2-1) = 0;
    TOK(s4(1,1)*2-1,s4(1,1)*2-1) = 1;
    TOF(s4(1,1)*2-1,1) = 0;
    for i = 2:size(s3,1)
        TOK(s3(i,1)*2,:) = 0;
        TOK(:,s3(i,1)*2) = 0;
        TOK(s3(i,1)*2,s3(i,1)*2) = 1;
        TOF(s3(i,1)*2,1) = 0;
    end
    for i = 2:size(s1,1)
        TOK(s1(i,1)*2,s1(i,1)*2) = C*TOK(s1(i,1)*2,s1(i,1)*2);
        TOF(s1(i,1)*2,1) = TOK(s1(i,1)*2,s1(i,1)*2)*((i-1)/(size(s1,1)-1));
    end
    for i = 2:size(s2,1)
        TOK(s2(i,1)*2,s2(i,1)*2) = C*TOK(s2(i,1)*2,s2(i,1)*2);
        TOF(s2(i,1)*2,1) = TOK(s2(i,1)*2,s2(i,1)*2)*(1-(i-1)/(size(s2,1)-1));
    end
    V = TOK \ TOF;
    ON2y = V(2:2:end);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % v3
    TOK = OK;
    TOF = OF;
    for i = 2:size(s1,1)
        TOK(s1(i,1)*2,:) = 0;
        TOK(:,s1(i,1)*2) = 0;
        TOK(s1(i,1)*2,s1(i,1)*2) = 1;
        TOF(s1(i,1)*2,1) = 0;
    end
    TOK(s1(1,1)*2-1,:) = 0;
    TOK(:,s1(1,1)*2-1) = 0;
    TOK(s1(1,1)*2-1,s1(1,1)*2-1) = 1;
    TOF(s1(1,1)*2-1,1) = 0;
    for i = 2:size(s4,1)
        TOK(s4(i,1)*2,:) = 0;
        TOK(:,s4(i,1)*2) = 0;
        TOK(s4(i,1)*2,s4(i,1)*2) = 1;
        TOF(s4(i,1)*2,1) = 0;
    end
    for i = 2:size(s2,1)
        TOK(s2(i,1)*2,s2(i,1)*2) = C*TOK(s2(i,1)*2,s2(i,1)*2);
        TOF(s2(i,1)*2,1) = TOK(s2(i,1)*2,s2(i,1)*2)*((i-1)/(size(s2,1)-1));
    end
    for i = 2:size(s3,1)
        TOK(s3(i,1)*2,s3(i,1)*2) = C*TOK(s3(i,1)*2,s3(i,1)*2);
        TOF(s3(i,1)*2,1) = TOK(s3(i,1)*2,s3(i,1)*2)*(1-(i-1)/(size(s3,1)-1));
    end
    V = TOK \ TOF;
    ON3y = V(2:2:end);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                                                                                                                                                                                                                                                       %%%%%%%%%%%%%%
    % v4
    TOK = OK;
    TOF = OF;
    for i = 2:size(s2,1)
        TOK(s2(i,1)*2,:) = 0;
        TOK(:,s2(i,1)*2) = 0;
        TOK(s2(i,1)*2,s2(i,1)*2) = 1;
        TOF(s2(i,1)*2,1) = 0;
    end
    TOK(s2(1,1)*2-1,:) = 0;
    TOK(:,s2(1,1)*2-1) = 0;
    TOK(s2(1,1)*2-1,s2(1,1)*2-1) = 1;
    TOF(s2(1,1)*2-1,1) = 0;
    for i = 2:size(s1,1)
        TOK(s1(i,1)*2,:) = 0;
        TOK(:,s1(i,1)*2) = 0;
        TOK(s1(i,1)*2,s1(i,1)*2) = 1;
        TOF(s1(i,1)*2,1) = 0;
    end
    for i = 2:size(s3,1)
        TOK(s3(i,1)*2,s3(i,1)*2) = C*TOK(s3(i,1)*2,s3(i,1)*2);
        TOF(s3(i,1)*2,1) = TOK(s3(i,1)*2,s3(i,1)*2)*((i-1)/(size(s3,1)-1));
    end
    for i = 2:size(s4,1)
        TOK(s4(i,1)*2,s4(i,1)*2) = C*TOK(s4(i,1)*2,s4(i,1)*2);
        TOF(s4(i,1)*2,1) = TOK(s4(i,1)*2,s4(i,1)*2)*(1-(i-1)/(size(s4,1)-1));
    end
    V = TOK \ TOF;
    ON4y = V(2:2:end);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    mc = zeros(4,4);
    for i = 1:4
        i1 = sub2os(VNos(i,1),1);
        mc(1,i) = ON1y(i1);
        mc(2,i) = ON2y(i1);
        mc(3,i) = ON3y(i1);
        mc(4,i) = ON4y(i1);
    end
    mc = inv(mc);
    SN1y = zeros(SNN,1);
    SN2y = zeros(SNN,1);
    SN3y = zeros(SNN,1);
    SN4y = zeros(SNN,1);
    for i = 1:SNN
        i1 = sub2os(i,1);
        SN1y(i,1) = mc(1,1)*ON1y(i1)+mc(1,2)*ON2y(i1)+mc(1,3)*ON3y(i1)+mc(1,4)*ON4y(i1);
        SN2y(i,1) = mc(2,1)*ON1y(i1)+mc(2,2)*ON2y(i1)+mc(2,3)*ON3y(i1)+mc(2,4)*ON4y(i1);
        SN3y(i,1) = mc(3,1)*ON1y(i1)+mc(3,2)*ON2y(i1)+mc(3,3)*ON3y(i1)+mc(3,4)*ON4y(i1);
        SN4y(i,1) = mc(4,1)*ON1y(i1)+mc(4,2)*ON2y(i1)+mc(4,3)*ON3y(i1)+mc(4,4)*ON4y(i1);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    % construction of the numerical shape function of the sub-grid mesh
    Xmin = min(S_Nodes(:,1)); 
    Ymin = min(S_Nodes(:,2));
    Xmax = max(S_Nodes(:,1));
    Ymax = max(S_Nodes(:,2));
    edge1 = find(abs(S_Nodes(:,2)-Ymin)<sy);
    edge2 = find(abs(S_Nodes(:,1)-Xmax)<sx);
    edge3 = find(abs(S_Nodes(:,2)-Ymax)<sy);
    edge4 = find(abs(S_Nodes(:,1)-Xmin)<sx);
    side = [edge1(:,1) S_Nodes(edge1(:,1),1)];
    side = sortrows(side,2);
    s1 = side(:,1);             
    side = [edge2(:,1) S_Nodes(edge2(:,1),2)];
    side = sortrows(side,2);
    s2 = side(:,1);            
    side = [edge3(:,1) S_Nodes(edge3(:,1),1)];
    side = sortrows(side,-2);
    s3 = side(:,1);             
    side = [edge4(:,1) S_Nodes(edge4(:,1),2)];
    side = sortrows(side,-2);
    s4 = side(:,1);              
    
    ns1 = size(s1,1);
    ns2 = size(s2,1);
    ns3 = size(s3,1);
    ns4 = size(s4,1);
    
    SF = sparse(2*SNN,1);
    SN = zeros(2*SNN,8);
    
    EK_Sub_Elems = zeros(64,SNE);
    for ise = 1:SNE
        E = SE(ise);
        D = Get_D(E,mu);
        EK = StiffnessMatrix_FineElement(S_Nodes(S_Elems(ise,:),:),D,t);
        EK_Sub_Elems(:,ise) = EK(:);
    end
    edofMat_sub = kron(S_Elems,[2,2])+repmat([-1,0],SNE,4);
    iK_sub = reshape(kron(edofMat_sub,ones(8,1))',64*SNE,1);
    jK_sub = reshape(kron(edofMat_sub,ones(1,8))',64*SNE,1);   
    sK_sub = reshape(EK_Sub_Elems,64*SNE,1);
    SK = sparse(iK_sub,jK_sub,sK_sub); 
    SK = (SK+SK')/2;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % su1
    TSK = SK;
    TSF = SF;   
    for i = 2:ns1
        nod = s1(i,1);
        dofx = nod*2-1;
        dofy = nod*2;
        val = SN1x(nod,1);
        TSK(dofx,dofx) = C*TSK(dofx,dofx);
        TSF(dofx,1) = TSK(dofx,dofx)*val;
        
        TSK(dofy,:) = 0;
        TSK(:,dofy) = 0;
        TSK(dofy,dofy) = 1;
        TSF(dofy,1) = 0;
    end
    for i = 2:ns2
        nod = s2(i,1);
        dofx = nod*2-1;
        dofy = nod*2;
        val = SN1x(nod,1);
        TSK(dofx,dofx) = C*TSK(dofx,dofx);
        TSF(dofx,1) = TSK(dofx,dofx)*val;
        
        TSK(dofy,:) = 0;
        TSK(:,dofy) = 0;
        TSK(dofy,dofy) = 1;
        TSF(dofy,1) = 0;    
    end    
    for i = 2:ns3
        nod = s3(i,1);
        dofx = nod*2-1;
        dofy = nod*2;
        val = SN1x(nod,1);
        TSK(dofx,dofx) = C*TSK(dofx,dofx);
        TSF(dofx,1) = TSK(dofx,dofx)*val;
        
        TSK(dofy,:) = 0;
        TSK(:,dofy) = 0;
        TSK(dofy,dofy) = 1;
        TSF(dofy,1) = 0;    
    end    
    for i = 2:ns4
        nod = s4(i,1);
        dofx = nod*2-1;
        dofy = nod*2;
        val = SN1x(nod,1);
        TSK(dofx,dofx) = C*TSK(dofx,dofx);
        TSF(dofx,1) = TSK(dofx,dofx)*val;
        
        TSK(dofy,:) = 0;
        TSK(:,dofy) = 0;
        TSK(dofy,dofy) = 1;
        TSF(dofy,1) = 0;    
    end    
    U = TSK \ TSF;
    SN(:,1) = U;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % su2
    TSK = SK;
    TSF = SF;   
    for i = 2:ns1
        nod = s1(i,1);
        dofx = nod*2-1;
        dofy = nod*2;
        val = SN2x(nod,1);
        TSK(dofx,dofx) = C*TSK(dofx,dofx);
        TSF(dofx,1) = TSK(dofx,dofx)*val;
        
        TSK(dofy,:) = 0;
        TSK(:,dofy) = 0;
        TSK(dofy,dofy) = 1;
        TSF(dofy,1) = 0;
    end
    for i = 2:ns2
        nod = s2(i,1);
        dofx = nod*2-1;
        dofy = nod*2;
        val = SN2x(nod,1);
        TSK(dofx,dofx) = C*TSK(dofx,dofx);
        TSF(dofx,1) = TSK(dofx,dofx)*val;
        
        TSK(dofy,:) = 0;
        TSK(:,dofy) = 0;
        TSK(dofy,dofy) = 1;
        TSF(dofy,1) = 0;   
    end    
    for i = 2:ns3
        nod = s3(i,1);
        dofx = nod*2-1;
        dofy = nod*2;
        val = SN2x(nod,1);
        TSK(dofx,dofx) = C*TSK(dofx,dofx);
        TSF(dofx,1) = TSK(dofx,dofx)*val;
        
        TSK(dofy,:) = 0;
        TSK(:,dofy) = 0;
        TSK(dofy,dofy) = 1;
        TSF(dofy,1) = 0;    
    end    
    for i = 2:ns4
        nod = s4(i,1);
        dofx = nod*2-1;
        dofy = nod*2;
        val = SN2x(nod,1);
        TSK(dofx,dofx) = C*TSK(dofx,dofx);
        TSF(dofx,1) = TSK(dofx,dofx)*val;
        
        TSK(dofy,:) = 0;
        TSK(:,dofy) = 0;
        TSK(dofy,dofy) = 1;
        TSF(dofy,1) = 0;
    end    
    U = TSK \ TSF;
    SN(:,3) = U;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % su3
    TSK = SK;
    TSF = SF;   
    for i = 2:ns1
        nod = s1(i,1);
        dofx = nod*2-1;
        dofy = nod*2;
        val = SN3x(nod,1);
        TSK(dofx,dofx) = C*TSK(dofx,dofx);
        TSF(dofx,1) = TSK(dofx,dofx)*val;
        
        TSK(dofy,:) = 0;
        TSK(:,dofy) = 0;
        TSK(dofy,dofy) = 1;
        TSF(dofy,1) = 0;
    end
    for i = 2:ns2
        nod = s2(i,1);
        dofx = nod*2-1;
        dofy = nod*2;
        val = SN3x(nod,1);
        TSK(dofx,dofx) = C*TSK(dofx,dofx);
        TSF(dofx,1) = TSK(dofx,dofx)*val;
        
        TSK(dofy,:) = 0;
        TSK(:,dofy) = 0;
        TSK(dofy,dofy) = 1;
        TSF(dofy,1) = 0;
    end    
    for i = 2:ns3
        nod = s3(i,1);
        dofx = nod*2-1;
        dofy = nod*2;
        val = SN3x(nod,1);
        TSK(dofx,dofx) = C*TSK(dofx,dofx);
        TSF(dofx,1) = TSK(dofx,dofx)*val;
        
        TSK(dofy,:) = 0;
        TSK(:,dofy) = 0;
        TSK(dofy,dofy) = 1;
        TSF(dofy,1) = 0;
    end    
    for i = 2:ns4
        nod = s4(i,1);
        dofx = nod*2-1;
        dofy = nod*2;
        val = SN3x(nod,1);
        TSK(dofx,dofx) = C*TSK(dofx,dofx);
        TSF(dofx,1) = TSK(dofx,dofx)*val;
        
        TSK(dofy,:) = 0;
        TSK(:,dofy) = 0;
        TSK(dofy,dofy) = 1;
        TSF(dofy,1) = 0;
    end    
    U = TSK \ TSF;
    SN(:,5) = U;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % su4
    TSK = SK;
    TSF = SF;   
    for i = 2:ns1
        nod = s1(i,1);
        dofx = nod*2-1;
        dofy = nod*2;
        val = SN4x(nod,1);
        TSK(dofx,dofx) = C*TSK(dofx,dofx);
        TSF(dofx,1) = TSK(dofx,dofx)*val;
        
        TSK(dofy,:) = 0;
        TSK(:,dofy) = 0;
        TSK(dofy,dofy) = 1;
        TSF(dofy,1) = 0;
    end
    for i = 2:ns2
        nod = s2(i,1);
        dofx = nod*2-1;
        dofy = nod*2;
        val = SN4x(nod,1);
        TSK(dofx,dofx) = C*TSK(dofx,dofx);
        TSF(dofx,1) = TSK(dofx,dofx)*val;
        
        TSK(dofy,:) = 0;
        TSK(:,dofy) = 0;
        TSK(dofy,dofy) = 1;
        TSF(dofy,1) = 0;
    end    
    for i = 2:ns3
        nod = s3(i,1);
        dofx = nod*2-1;
        dofy = nod*2;
        val = SN4x(nod,1);
        TSK(dofx,dofx) = C*TSK(dofx,dofx);
        TSF(dofx,1) = TSK(dofx,dofx)*val;
        
        TSK(dofy,:) = 0;
        TSK(:,dofy) = 0;
        TSK(dofy,dofy) = 1;
        TSF(dofy,1) = 0;
    end    
    for i = 2:ns4
        nod = s4(i,1);
        dofx = nod*2-1;
        dofy = nod*2;
        val = SN4x(nod,1);
        TSK(dofx,dofx) = C*TSK(dofx,dofx);
        TSF(dofx,1) = TSK(dofx,dofx)*val;
        
        TSK(dofy,:) = 0;
        TSK(:,dofy) = 0;
        TSK(dofy,dofy) = 1;
        TSF(dofy,1) = 0;
    end    
    U = TSK \ TSF;
    SN(:,7) = U;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % sv1
    TSK = SK;
    TSF = SF;   
    for i = 2:ns1
        nod = s1(i,1);
        dofx = nod*2-1;
        dofy = nod*2;
        val = SN1y(nod,1);
        TSK(dofy,dofy) = C*TSK(dofy,dofy);
        TSF(dofy,1) = TSK(dofy,dofy)*val;
        
        TSK(dofx,:) = 0;
        TSK(:,dofx) = 0;
        TSK(dofx,dofx) = 1;
        TSF(dofx,1) = 0;
    end
    for i = 2:ns2
        nod = s2(i,1);
        dofx = nod*2-1;
        dofy = nod*2;
        val = SN1y(nod,1);
        TSK(dofy,dofy) = C*TSK(dofy,dofy);
        TSF(dofy,1) = TSK(dofy,dofy)*val;
        
        TSK(dofx,:) = 0;
        TSK(:,dofx) = 0;
        TSK(dofx,dofx) = 1;
        TSF(dofx,1) = 0;    
    end    
    for i = 2:ns3
        nod = s3(i,1);
        dofx = nod*2-1;
        dofy = nod*2;
        val = SN1y(nod,1);
        TSK(dofy,dofy) = C*TSK(dofy,dofy);
        TSF(dofy,1) = TSK(dofy,dofy)*val;
        
        TSK(dofx,:) = 0;
        TSK(:,dofx) = 0;
        TSK(dofx,dofx) = 1;
        TSF(dofx,1) = 0;      
    end    
    for i = 2:ns4
        nod = s4(i,1);
        dofx = nod*2-1;
        dofy = nod*2;
        val = SN1y(nod,1);
        TSK(dofy,dofy) = C*TSK(dofy,dofy);
        TSF(dofy,1) = TSK(dofy,dofy)*val;
        
        TSK(dofx,:) = 0;
        TSK(:,dofx) = 0;
        TSK(dofx,dofx) = 1;
        TSF(dofx,1) = 0;  
    end    
    U = TSK \ TSF;
    SN(:,2) = U;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % sv2
    TSK = SK;
    TSF = SF;   
    for i = 2:ns1
        nod = s1(i,1);
        dofx = nod*2-1;
        dofy = nod*2;
        val = SN2y(nod,1);
        TSK(dofy,dofy) = C*TSK(dofy,dofy);
        TSF(dofy,1) = TSK(dofy,dofy)*val;
        
        TSK(dofx,:) = 0;
        TSK(:,dofx) = 0;
        TSK(dofx,dofx) = 1;
        TSF(dofx,1) = 0;
    end
    for i = 2:ns2
        nod = s2(i,1);
        dofx = nod*2-1;
        dofy = nod*2;
        val = SN2y(nod,1);
        TSK(dofy,dofy) = C*TSK(dofy,dofy);
        TSF(dofy,1) = TSK(dofy,dofy)*val;
        
        TSK(dofx,:) = 0;
        TSK(:,dofx) = 0;
        TSK(dofx,dofx) = 1;
        TSF(dofx,1) = 0;    
    end    
    for i = 2:ns3
        nod = s3(i,1);
        dofx = nod*2-1;
        dofy = nod*2;
        val = SN2y(nod,1);
        TSK(dofy,dofy) = C*TSK(dofy,dofy);
        TSF(dofy,1) = TSK(dofy,dofy)*val;
        
        TSK(dofx,:) = 0;
        TSK(:,dofx) = 0;
        TSK(dofx,dofx) = 1;
        TSF(dofx,1) = 0;       
    end    
    for i = 2:ns4
        nod = s4(i,1);
        dofx = nod*2-1;
        dofy = nod*2;
        val = SN2y(nod,1);
        TSK(dofy,dofy) = C*TSK(dofy,dofy);
        TSF(dofy,1) = TSK(dofy,dofy)*val;
        
        TSK(dofx,:) = 0;
        TSK(:,dofx) = 0;
        TSK(dofx,dofx) = 1;
        TSF(dofx,1) = 0;      
    end    
    U = TSK \ TSF;
    SN(:,4) = U;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % sv3
    TSK = SK;
    TSF = SF;   
    for i = 2:ns1
        nod = s1(i,1);
        dofx = nod*2-1;
        dofy = nod*2;
        val = SN3y(nod,1);
        TSK(dofy,dofy) = C*TSK(dofy,dofy);
        TSF(dofy,1) = TSK(dofy,dofy)*val;
        
        TSK(dofx,:) = 0;
        TSK(:,dofx) = 0;
        TSK(dofx,dofx) = 1;
        TSF(dofx,1) = 0;
    end
    for i = 2:ns2
        nod = s2(i,1);
        dofx = nod*2-1;
        dofy = nod*2;
        val = SN3y(nod,1);
        TSK(dofy,dofy) = C*TSK(dofy,dofy);
        TSF(dofy,1) = TSK(dofy,dofy)*val;
        
        TSK(dofx,:) = 0;
        TSK(:,dofx) = 0;
        TSK(dofx,dofx) = 1;
        TSF(dofx,1) = 0;   
    end    
    for i = 2:ns3
        nod = s3(i,1);
        dofx = nod*2-1;
        dofy = nod*2;
        val = SN3y(nod,1);
        TSK(dofy,dofy) = C*TSK(dofy,dofy);
        TSF(dofy,1) = TSK(dofy,dofy)*val;
        
        TSK(dofx,:) = 0;
        TSK(:,dofx) = 0;
        TSK(dofx,dofx) = 1;
        TSF(dofx,1) = 0;      
    end    
    for i = 2:ns4
        nod = s4(i,1);
        dofx = nod*2-1;
        dofy = nod*2;
        val = SN3y(nod,1);
        TSK(dofy,dofy) = C*TSK(dofy,dofy);
        TSF(dofy,1) = TSK(dofy,dofy)*val;
        
        TSK(dofx,:) = 0;
        TSK(:,dofx) = 0;
        TSK(dofx,dofx) = 1;
        TSF(dofx,1) = 0;       
    end    
    U = TSK \ TSF;
    SN(:,6) = U;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % sv4
    TSK = SK;
    TSF = SF;   
    for i = 2:ns1
        nod = s1(i,1);
        dofx = nod*2-1;
        dofy = nod*2;
        val = SN4y(nod,1);
        TSK(dofy,dofy) = C*TSK(dofy,dofy);
        TSF(dofy,1) = TSK(dofy,dofy)*val;
        
        TSK(dofx,:) = 0;
        TSK(:,dofx) = 0;
        TSK(dofx,dofx) = 1;
        TSF(dofx,1) = 0;
    end
    for i = 2:ns2
        nod = s2(i,1);
        dofx = nod*2-1;
        dofy = nod*2;
        val = SN4y(nod,1);
        TSK(dofy,dofy) = C*TSK(dofy,dofy);
        TSF(dofy,1) = TSK(dofy,dofy)*val;
        
        TSK(dofx,:) = 0;
        TSK(:,dofx) = 0;
        TSK(dofx,dofx) = 1;
        TSF(dofx,1) = 0; 
    end    
    for i = 2:ns3
        nod = s3(i,1);
        dofx = nod*2-1;
        dofy = nod*2;
        val = SN4y(nod,1);
        TSK(dofy,dofy) = C*TSK(dofy,dofy);
        TSF(dofy,1) = TSK(dofy,dofy)*val;
        
        TSK(dofx,:) = 0;
        TSK(:,dofx) = 0;
        TSK(dofx,dofx) = 1;
        TSF(dofx,1) = 0;      
    end    
    for i = 2:ns4
        nod = s4(i,1);
        dofx = nod*2-1;
        dofy = nod*2;
        val = SN4y(nod,1);
        TSK(dofy,dofy) = C*TSK(dofy,dofy);
        TSF(dofy,1) = TSK(dofy,dofy)*val;
        
        TSK(dofx,:) = 0;
        TSK(:,dofx) = 0;
        TSK(dofx,dofx) = 1;
        TSF(dofx,1) = 0;      
    end    
    U = TSK \ TSF;
    SN(:,8) = U;

    %% Calculate the macroscopic equivalent matrix
    % 
    Q = sparse(2*SNN,1);
    
    % external force is imposed on the right edge
    q = -1/4;
    % X direction
%     Q((2*s2(1)-1),1) = q/(ns2-1)/2;
%     Q((2*s2(ns2)-1),1) = q/(ns2-1)/2;
%     Q((2*s2(2:(ns2-1))-1),1) = q/(ns2-1);
    % Y direction
    Q((2*s2(1)),1) = q/(ns2-1)/2;
    Q((2*s2(ns2)),1) = q/(ns2-1)/2;
    Q((2*s2(2:(ns2-1))),1) = q/(ns2-1);


    % Y direction, external force is imposed on the top edge
%     q = -1/2;
%     Q((2*s3(1)),1) = q/(ns3-1)/2;
%     Q((2*s3(ns3)),1) = q/(ns3-1)/2;
%     Q((2*s3(2:(ns3-1))),1) = q/(ns3-1);

    % X direction, external force is imposed on the left edge
%     q = 1/18;
%     Q((2*s4(1)-1),1) = q/(ns4-1)/2;
%     Q((2*s4(ns4)-1),1) = q/(ns4-1)/2;
%     Q((2*s4(2:(ns4-1))-1),1) = q/(ns4-1);
    
    
    CEKs = SN'*SK*SN;
    CEFs = SN'*Q;
    
    
        
end