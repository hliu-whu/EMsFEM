function VNos = VerticesNodesNo(TO_Nodes,Xmin,Ymin,Xmax,Ymax,sx,sy)

% Get the node numbers of the four corners of the oversampling mesh
    VNos = zeros(4,1);
    VNos(1,1) = find((abs(TO_Nodes(:,1)-Xmin)<sx) & (abs(TO_Nodes(:,2)-Ymin)<sy));
    VNos(2,1) = find((abs(TO_Nodes(:,1)-Xmax)<sx) & (abs(TO_Nodes(:,2)-Ymin)<sy));
    VNos(3,1) = find((abs(TO_Nodes(:,1)-Xmax)<sx) & (abs(TO_Nodes(:,2)-Ymax)<sy));
    VNos(4,1) = find((abs(TO_Nodes(:,1)-Xmin)<sx) & (abs(TO_Nodes(:,2)-Ymax)<sy));

end