function phi_fine = GetFineRes_Hetero(C_Nodes,C_Elems,S_Nodes,phi_coarse,SN_Cells)
% for non-periodic structures

    mx = 1;
    my = 1;
    sx = 1e-6;
    sy = 1e-6;
    F_Nodes = load('..\data\Fine_Nodes.dat');

    FNN = size(F_Nodes,1);
%     CNN = size(C_Nodes,1);
    CNE = size(C_Elems,1);
    SNN = size(S_Nodes,1);
    TS_Nodes = zeros(SNN,2);
    cn = zeros(FNN,1);
    N = size(phi_coarse,2);
    phi_fine = zeros(2*FNN,N);
    Cndofs = zeros(1,4*(mx+my));
    for ice = 1:CNE
        CEN1 = C_Nodes(C_Elems(ice,1),:);
        TS_Nodes(:,1) = S_Nodes(:,1) + CEN1(1,1);
        TS_Nodes(:,2) = S_Nodes(:,2) + CEN1(1,2);
        Xmax = max(TS_Nodes(:,1));
        Xmin = min(TS_Nodes(:,1));
        Ymax = max(TS_Nodes(:,2));
        Ymin = min(TS_Nodes(:,2));
        fdofs = zeros(1,2*SNN);
        SN = SN_Cells{ice,1}; % load(strcat('..\data\SN',int2str(ice),'.dat'));
        for isn = 1:SNN
            nx = TS_Nodes(isn,1); 
            ny = TS_Nodes(isn,2);
            if((abs(ny-Ymin)<sy) || (abs(ny-Ymax)<sy) || (abs(nx-Xmin)<sx) || (abs(nx-Xmax)<sx))
                marker = 1;
            else
                marker = 0;
            end
            node = find((abs(nx-F_Nodes(:,1))<sx) & (abs(ny-F_Nodes(:,2))<sy));
            [a,b] = size(node);
            if(a*b~=1)
                disp('Error in Function GetFineModes!!!');
            else
                fdofs(1,2*isn-1) = 2*node-1;
                fdofs(1,2*isn)   = 2*node;
                if(marker)
                    cn(node,1) = cn(node,1) + 1;
                end
            end
        end
        Cndofs(1,1:2:end) = C_Elems(ice,:)*2-1;         Cndofs(1,2:2:end) = C_Elems(ice,:)*2;
        Un = phi_coarse(Cndofs,:);

        phi_fine(fdofs,:) = phi_fine(fdofs,:) + SN*Un;

    end
    
    for ifn = 1:FNN
        if(cn(ifn,1)>1)
            phi_fine(2*ifn-1,:) = phi_fine(2*ifn-1,:)/cn(ifn,1);
            phi_fine(2*ifn,:)   = phi_fine(2*ifn,:)/cn(ifn,1);
        end
    end

end