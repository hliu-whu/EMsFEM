function [MatMappingOS,Sub2OSElems] = MatMappingOS_fun(C_Nodes,C_Elems,sxyz)
% Material mapping among the fine-scale mesh, the sub-grid mesh and the oversampling mesh

    F_Nodes = load('../data/Fine_Nodes.dat');
    F_Elems = load('../data/Fine_Elements.dat');
    F_Mods  = load('../data/Fine_Modulus.dat');
    O_Nodes = load('../data/OS_Nodes.dat');
    O_Elems = load('../data/OS_Elements.dat');
    S_Nodes = load('../data/Sub_Nodes.dat');
    S_Elems = load('../data/Sub_Elements.dat');
    CNE = size(C_Elems,1);
    ONE = size(O_Elems,1);
    SNE = size(S_Elems,1);
    FNE = size(F_Elems,1);
    FXmin = min(F_Nodes(:,1));
    FXmax = max(F_Nodes(:,1));
    FYmin = min(F_Nodes(:,2));
    FYmax = max(F_Nodes(:,2));

    MatMappingOS = zeros(ONE,CNE);
    Sub2OSElems = zeros(SNE,1);
    CenterOfFineElems = zeros(FNE,2); % The coordinates of the center point of elements in the fine-scale mesh in the x and y direction, respectively
    CenterOfFineElems(:,1) = ( F_Nodes(F_Elems(:,1),1) + F_Nodes(F_Elems(:,2),1) + F_Nodes(F_Elems(:,3),1) + F_Nodes(F_Elems(:,4),1) )/4;
    CenterOfFineElems(:,2) = ( F_Nodes(F_Elems(:,1),2) + F_Nodes(F_Elems(:,2),2) + F_Nodes(F_Elems(:,3),2) + F_Nodes(F_Elems(:,4),2) )/4;

    CenterOfOSElems = zeros(ONE,2); % The coordinates of the center point of elements in the oversampling mesh in the x and y direction, respectively
    CenterOfOSElems(:,1) = ( O_Nodes(O_Elems(:,1),1) + O_Nodes(O_Elems(:,2),1) + O_Nodes(O_Elems(:,3),1) + O_Nodes(O_Elems(:,4),1) )/4;
    CenterOfOSElems(:,2) = ( O_Nodes(O_Elems(:,1),2) + O_Nodes(O_Elems(:,2),2) + O_Nodes(O_Elems(:,3),2) + O_Nodes(O_Elems(:,4),2) )/4;
    
    CenterOfSubElems = zeros(SNE,2); % The coordinates of the center point of elements in the sub-grid mesh in the x and y direction, respectively
    CenterOfSubElems(:,1) = ( S_Nodes(S_Elems(:,1),1) + S_Nodes(S_Elems(:,2),1) + S_Nodes(S_Elems(:,3),1) + S_Nodes(S_Elems(:,4),1) )/4;
    CenterOfSubElems(:,2) = ( S_Nodes(S_Elems(:,1),2) + S_Nodes(S_Elems(:,2),2) + S_Nodes(S_Elems(:,3),2) + S_Nodes(S_Elems(:,4),2) )/4;
    
    TO_Center = zeros(ONE,2);
    Xmax = max(S_Nodes(:,1));
    Ymax = max(S_Nodes(:,2));
    
    for ice = 1:CNE
        CN1 = C_Nodes(C_Elems(ice,1),:); % Coordinates of the first node of the coarse element
        TO_Center(:,1) = CenterOfOSElems(:,1) + CN1(1,1) - Xmax;
        TO_Center(:,2) = CenterOfOSElems(:,2) + CN1(1,2) - Ymax;
        CenterOSx = sum(TO_Center(:,1))/ONE;
        
        if( (CenterOSx>FXmin) && (CenterOSx<FXmax) )
            id = TO_Center(:,1) < FXmin;
            if(length(id)>=1)
                TO_Center(id,1) = 2*FXmin - TO_Center(id,1);
            end
            id = TO_Center(:,1) > FXmax;
            if(length(id)>=1)
                TO_Center(id,1) = 2*FXmax - TO_Center(id,1);
            end
            id = TO_Center(:,2) < FYmin;
            if(length(id)>=1)
                TO_Center(id,2) = 2*FYmin - TO_Center(id,2);
            end
            id = TO_Center(:,2) > FYmax;
            if(length(id)>=1)
                TO_Center(id,2) = 2*FYmax - TO_Center(id,2);
            end
        else
            error('CenterOSx!!!');
        end

        for ioe = 1:ONE
            id = (abs(TO_Center(ioe,1)-CenterOfFineElems(:,1))<sxyz) & (abs(TO_Center(ioe,2)-CenterOfFineElems(:,2))<sxyz);
            ife = find(id);
            if(length(ife)==1)
                MatMappingOS(ioe,ice) = F_Mods(ife,1);
            else
                error('yi wai MatMappingOS!!!');
            end
        end
    end
    
    TS_Center = zeros(SNE,1);
    TS_Center(:,1) = CenterOfSubElems(:,1) + Xmax;
    TS_Center(:,2) = CenterOfSubElems(:,2) + Ymax;
    
    for ise = 1:SNE
        id = (abs(TS_Center(ise,1)-CenterOfOSElems(:,1))<sxyz) & (abs(TS_Center(ise,2)-CenterOfOSElems(:,2))<sxyz);
        ioe = find(id);
        if(length(ioe)==1)
            Sub2OSElems(ise,1) = ioe;
        else
            error('yi wai Sub2OSElems!!!');
        end
    end
    
    
    
end