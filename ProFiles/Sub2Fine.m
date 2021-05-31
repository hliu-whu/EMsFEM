function [map_sub2fine] = Sub2Fine(ice)

% Generate element mapping from the sub-grid mesh to the fine-scale mesh
% This specifies that the lower left corner node of the coarse element is the first node

    C_Elems = load('..\data\Coarse_Elements.dat');
    C_Nodes = load('..\data\Coarse_Nodes.dat');
    S_Nodes = load('..\data\Sub_Nodes.dat');
    F_Nodes = load('..\data\Fine_Nodes.dat');
    
    sxyz = 1e-6;

    SNN = size(S_Nodes,1);

    map_sub2fine = zeros(SNN,1);
    
    hx = F_Nodes(:,1);
    hy = F_Nodes(:,2);
    
    
        TS_Nodes = S_Nodes;
        TS_Nodes(:,1) = TS_Nodes(:,1) + C_Nodes(C_Elems(ice,1),1);
        TS_Nodes(:,2) = TS_Nodes(:,2) + C_Nodes(C_Elems(ice,1),2);
        for isn = 1:SNN
            sx = TS_Nodes(isn,1);
            sy = TS_Nodes(isn,2);
            nx = (abs(hx-sx)<sxyz) & (abs(hy-sy)<sxyz);
            mx = find(nx);
            [a,b] = size(mx);
            if(a*b~=1)
                disp('Error1 in the function Sub2Fine!');
            else
                map_sub2fine(isn,1) = mx;
            end
        end
        clear TS_Nodes;
        
%         TO_Nodes = O_Nodes;
%         TO_Nodes(:,1) = TO_Nodes(:,1) + C_Nodes(C_Elems(ice,1),1) + shiftxy(elementtype(ice,1),1);
%         TO_Nodes(:,2) = TO_Nodes(:,2) + C_Nodes(C_Elems(ice,1),2) + shiftxy(elementtype(ice,1),2);
%         for ioe = 1:ONE
%             sx = sum(TO_Nodes(O_Elems(ioe,:),1))/4;
%             sy = sum(TO_Nodes(O_Elems(ioe,:),2))/4;
%             nx = (abs(hx-sx)<scalex) & (abs(hy-sy)<scaley);
%             mx = find(nx);
%             [a,b] = size(mx);
%             if(a*b~=1)
%                 disp('Error2 in the function SubOs2Fine!');
%             else
%                 map_os2fine(ioe,ice) = mx;
%             end            
%         end
        
%         if(mod(ice,100)==0)
%             fprintf('ice = %6d\n',ice);
%         end


end