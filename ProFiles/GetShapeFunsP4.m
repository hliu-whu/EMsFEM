function [SN,CEKs,CEFs] = GetShapeFunsP4(S_Nodes,S_Elems,SK,sx,sy)

    C = abs(max(SK(:)))*1.0e6;

    SNN = size(S_Nodes,1);
    Xmin = min(S_Nodes(:,1)); % four corners
    Ymin = min(S_Nodes(:,2));
    Xmax = max(S_Nodes(:,1));
    Ymax = max(S_Nodes(:,2));
    edge1 = find(abs(S_Nodes(:,2)-Ymin)<sy);
    edge2 = find(abs(S_Nodes(:,1)-Xmax)<sx);
    edge3 = find(abs(S_Nodes(:,2)-Ymax)<sy);
    edge4 = find(abs(S_Nodes(:,1)-Xmin)<sx);
    side = [edge1(:,1) S_Nodes(edge1(:,1),1)];
    side = sortrows(side,2);
    s1 = side(:,1);              % nodes on the lower edge, ordering from left to right
    side = [edge2(:,1) S_Nodes(edge2(:,1),2)];
    side = sortrows(side,2);
    s2 = side(:,1);              % nodes on the right edge, ordering from bottom to top
    side = [edge3(:,1) S_Nodes(edge3(:,1),1)];
    side = sortrows(side,2);
    s3 = side(:,1);              % nodes on the upper edge, ordering from left to right
    side = [edge4(:,1) S_Nodes(edge4(:,1),2)];
    side = sortrows(side,2);
    s4 = side(:,1);              % nodes on the left edge, ordering from bottom to top
    
    ns1 = size(s1,1);
    ns2 = size(s2,1);
    ns3 = size(s3,1);
    ns4 = size(s4,1);

    SF = sparse(2*SNN,1);
    SN = zeros(2*SNN,8);

    % Construct the numerical shape functions by using the periodic boundary conditions
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % for the first macro node of coarse element in the x direction
    % u1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    TOK = SK;
    TOF = SF;
    bt1 = 1;
    bt2 = -1;
    Cbt = C*[bt1*bt1 bt1*bt2; bt2*bt1 bt2*bt2];

    TOK(s2(ns2,1)*2-1,s2(ns2,1)*2-1) = TOK(s2(ns2,1)*2-1,s2(ns2,1)*2-1) + C;
    TOF(s2(ns2,1)*2-1,1) = 0;
    TOK(s2(ns2,1)*2,s2(ns2,1)*2) = TOK(s2(ns2,1)*2,s2(ns2,1)*2) + C;
    TOF(s2(ns2,1)*2,1) = 0;
    
    S1x = linspace(1,0,ns1);
    for i = 1:ns1
        Vdofs = [s1(i,1)*2 s3(i,1)*2];
        TOK(Vdofs,Vdofs) = TOK(Vdofs,Vdofs) + Cbt;
        Udofs = [s1(i,1)*2-1 s3(i,1)*2-1];
        TOK(Udofs,Udofs) = TOK(Udofs,Udofs) + Cbt;

        bt0 = S1x(i);
        FCbt = C*[bt0*bt1;bt0*bt2];
        TOF(Udofs,1) = TOF(Udofs,1) + FCbt;
    end

    S4x = linspace(1,0,ns4);
    for i = 1:(ns4-1)
        Vdofs = [s4(i,1)*2 s2(i,1)*2];
        TOK(Vdofs,Vdofs) = TOK(Vdofs,Vdofs) + Cbt;
        Udofs = [s4(i,1)*2-1 s2(i,1)*2-1];
        TOK(Udofs,Udofs) = TOK(Udofs,Udofs) + Cbt;
        
        bt0 = S4x(i);
        FCbt = C*[bt0*bt1;bt0*bt2];
        TOF(Udofs,1) = TOF(Udofs,1) + FCbt;
    end
    U = TOK \ TOF;
    SN(:,1) = U(:,1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % for the first macro node of coarse element in the x direction
    % u2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    TOK = SK;
    TOF = SF;
    
    TOK(s3(ns3,1)*2-1,s3(ns3,1)*2-1) = C*TOK(s3(ns3,1)*2-1,s3(ns3,1)*2-1);
    TOF(s3(ns3,1)*2-1,1) = 0;
    TOK(s3(ns3,1)*2,s3(ns3,1)*2) = C*TOK(s3(ns3,1)*2,s3(ns3,1)*2);
    TOF(s3(ns3,1)*2,1) = 0;
    
    S1x = linspace(0,1,ns1);
    for i = 1:ns1
        Vdofs = [s1(i,1)*2 s3(i,1)*2];
        TOK(Vdofs,Vdofs) = TOK(Vdofs,Vdofs) + Cbt;
        Udofs = [s1(i,1)*2-1 s3(i,1)*2-1];
        TOK(Udofs,Udofs) = TOK(Udofs,Udofs) + Cbt;
        
        bt0 = S1x(i);
        FCbt = C*[bt0*bt1;bt0*bt2];
        TOF(Udofs,1) = TOF(Udofs,1) + FCbt;
    end
    
    S2x = linspace(1,0,ns2);
    for i = 2:ns2
        Vdofs = [s2(i,1)*2 s4(i,1)*2];
        TOK(Vdofs,Vdofs) = TOK(Vdofs,Vdofs) + Cbt;
        Udofs = [s2(i,1)*2-1 s4(i,1)*2-1];
        TOK(Udofs,Udofs) = TOK(Udofs,Udofs) + Cbt;
        
        bt0 = S2x(i);
        FCbt = C*[bt0*bt1;bt0*bt2];
        TOF(Udofs,1) = TOF(Udofs,1) + FCbt;
    end
    U = TOK \ TOF;
    SN(:,3) = U(:,1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % for the first macro node of coarse element in the x direction
    % u3
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    TOK = SK;
    TOF = SF;
    
    TOK(s1(1,1)*2-1,s1(1,1)*2-1) = C*TOK(s1(1,1)*2-1,s1(1,1)*2-1); 
    TOF(s1(1,1)*2-1,1) = 0;
    TOK(s1(1,1)*2,s1(1,1)*2) = C*TOK(s1(1,1)*2,s1(1,1)*2);
    TOF(s1(1,1)*2,1) = 0;
    
    S3x = linspace(0,1,ns3);
    for i = 1:ns3
        Vdofs = [s3(i,1)*2 s1(i,1)*2];
        TOK(Vdofs,Vdofs) = TOK(Vdofs,Vdofs) + Cbt;
        Udofs = [s3(i,1)*2-1 s1(i,1)*2-1];
        TOK(Udofs,Udofs) = TOK(Udofs,Udofs) + Cbt;
        
        bt0 = S3x(i);
        FCbt = C*[bt0*bt1;bt0*bt2];
        TOF(Udofs,1) = TOF(Udofs,1) + FCbt;
    end
    
    S2x = linspace(0,1,ns2);
    for i = 1:(ns2-1)
        Vdofs = [s2(i,1)*2 s4(i,1)*2];
        TOK(Vdofs,Vdofs) = TOK(Vdofs,Vdofs) + Cbt;
        Udofs = [s2(i,1)*2-1 s4(i,1)*2-1];
        TOK(Udofs,Udofs) = TOK(Udofs,Udofs) + Cbt;
        
        bt0 = S2x(i);
        FCbt = C*[bt0*bt1;bt0*bt2];
        TOF(Udofs,1) = TOF(Udofs,1) + FCbt;
    end
    U = TOK \ TOF;
    SN(:,5) = U(:,1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % for the first macro node of coarse element in the x direction
    % u4
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    TOK = SK;
    TOF = SF;
    
    TOK(s1(ns1,1)*2-1,s1(ns1,1)*2-1) = C*TOK(s1(ns1,1)*2-1,s1(ns1,1)*2-1); 
    TOF(s1(ns1,1)*2-1,1) = 0;
    TOK(s1(ns1,1)*2,s1(ns1,1)*2) = C*TOK(s1(ns1,1)*2,s1(ns1,1)*2);
    TOF(s1(ns1,1)*2,1) = 0;
    
    S3x = linspace(1,0,ns3);
    for i = 1:ns3
        Vdofs = [s3(i,1)*2 s1(i,1)*2];
        TOK(Vdofs,Vdofs) = TOK(Vdofs,Vdofs) + Cbt;
        Udofs = [s3(i,1)*2-1 s1(i,1)*2-1];
        TOK(Udofs,Udofs) = TOK(Udofs,Udofs) + Cbt;
        
        bt0 = S3x(i);
        FCbt = C*[bt0*bt1;bt0*bt2];
        TOF(Udofs,1) = TOF(Udofs,1) + FCbt;
    end
    
    S4x = linspace(0,1,ns4);
    for i = 1:(ns4-1)
        Vdofs = [s4(i,1)*2 s2(i,1)*2];
        TOK(Vdofs,Vdofs) = TOK(Vdofs,Vdofs) + Cbt;
        Udofs = [s4(i,1)*2-1 s2(i,1)*2-1];
        TOK(Udofs,Udofs) = TOK(Udofs,Udofs) + Cbt;
        
        bt0 = S4x(i);
        FCbt = C*[bt0*bt1;bt0*bt2];
        TOF(Udofs,1) = TOF(Udofs,1) + FCbt;
    end
    U = TOK \ TOF;
    SN(:,7) = U(:,1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % for the first macro node of coarse element in the x direction
    % v1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    TOK = SK;
    TOF = SF;
    
    TOK(s2(ns2,1)*2-1,s2(ns2,1)*2-1) = C*TOK(s2(ns2,1)*2-1,s2(ns2,1)*2-1); 
    TOF(s2(ns2,1)*2-1,1) = 0;
    TOK(s2(ns2,1)*2,s2(ns2,1)*2) = C*TOK(s2(ns2,1)*2,s2(ns2,1)*2);
    TOF(s2(ns2,1)*2,1) = 0;
    
    S1y = linspace(1,0,ns1);
    for i = 1:ns1
        Vdofs = [s1(i,1)*2 s3(i,1)*2];
        TOK(Vdofs,Vdofs) = TOK(Vdofs,Vdofs) + Cbt;
        Udofs = [s1(i,1)*2-1 s3(i,1)*2-1];
        TOK(Udofs,Udofs) = TOK(Udofs,Udofs) + Cbt;
        
        bt0 = S1y(i);
        FCbt = C*[bt0*bt1;bt0*bt2];
        TOF(Vdofs,1) = TOF(Vdofs,1) + FCbt;
    end
    
    S4y = linspace(1,0,ns4);
    for i = 2:ns4
        Vdofs = [s4(i,1)*2 s2(i,1)*2];
        TOK(Vdofs,Vdofs) = TOK(Vdofs,Vdofs) + Cbt;
        Udofs = [s4(i,1)*2-1 s2(i,1)*2-1];
        TOK(Udofs,Udofs) = TOK(Udofs,Udofs) + Cbt;
        
        bt0 = S4y(i);
        FCbt = C*[bt0*bt1;bt0*bt2];
        TOF(Vdofs,1) = TOF(Vdofs,1) + FCbt;
    end
    U = TOK \ TOF;
    SN(:,2) = U(:,1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % for the first macro node of coarse element in the x direction
    % v2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    TOK = SK;
    TOF = SF;
    
    TOK(s3(ns3,1)*2-1,s3(ns3,1)*2-1) = C*TOK(s3(ns3,1)*2-1,s3(ns3,1)*2-1); 
    TOF(s3(ns3,1)*2-1,1) = 0;
    TOK(s3(ns3,1)*2,s3(ns3,1)*2) = C*TOK(s3(ns3,1)*2,s3(ns3,1)*2);
    TOF(s3(ns3,1)*2,1) = 0;
    
    S1y = linspace(0,1,ns1);
    for i = 1:ns1
        Vdofs = [s1(i,1)*2 s3(i,1)*2];
        TOK(Vdofs,Vdofs) = TOK(Vdofs,Vdofs) + Cbt;
        Udofs = [s1(i,1)*2-1 s3(i,1)*2-1];
        TOK(Udofs,Udofs) = TOK(Udofs,Udofs) + Cbt;
        
        bt0 = S1y(i);
        FCbt = C*[bt0*bt1;bt0*bt2];
        TOF(Vdofs,1) = TOF(Vdofs,1) + FCbt;
    end
    
    S2y = linspace(1,0,ns2);
    for i = 2:ns2
        Vdofs = [s2(i,1)*2 s4(i,1)*2];
        TOK(Vdofs,Vdofs) = TOK(Vdofs,Vdofs) + Cbt;
        Udofs = [s2(i,1)*2-1 s4(i,1)*2-1];
        TOK(Udofs,Udofs) = TOK(Udofs,Udofs) + Cbt;
        
        bt0 = S2y(i);
        FCbt = C*[bt0*bt1;bt0*bt2];
        TOF(Vdofs,1) = TOF(Vdofs,1) + FCbt;
    end
    U = TOK \ TOF;
    SN(:,4) = U(:,1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % for the first macro node of coarse element in the x direction
    % v3
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    TOK = SK;
    TOF = SF;
    
    TOK(s1(1,1)*2-1,s1(1,1)*2-1) = C*TOK(s1(1,1)*2-1,s1(1,1)*2-1); 
    TOF(s1(1,1)*2-1,1) = 0;
    TOK(s1(1,1)*2,s1(1,1)*2) = C*TOK(s1(1,1)*2,s1(1,1)*2);
    TOF(s1(1,1)*2,1) = 0;
    
    S3y = linspace(0,1,ns3);
    for i = 1:ns3
        Vdofs = [s3(i,1)*2 s1(i,1)*2];
        TOK(Vdofs,Vdofs) = TOK(Vdofs,Vdofs) + Cbt;
        Udofs = [s3(i,1)*2-1 s1(i,1)*2-1];
        TOK(Udofs,Udofs) = TOK(Udofs,Udofs) + Cbt;
        
        bt0 = S3y(i);
        FCbt = C*[bt0*bt1;bt0*bt2];
        TOF(Vdofs,1) = TOF(Vdofs,1) + FCbt;
    end
    
    S2y = linspace(0,1,ns2);
    for i = 1:(ns2-1)
        Vdofs = [s2(i,1)*2 s4(i,1)*2];
        TOK(Vdofs,Vdofs) = TOK(Vdofs,Vdofs) + Cbt;
        Udofs = [s2(i,1)*2-1 s4(i,1)*2-1];
        TOK(Udofs,Udofs) = TOK(Udofs,Udofs) + Cbt;
        
        bt0 = S2y(i);
        FCbt = C*[bt0*bt1;bt0*bt2];
        TOF(Vdofs,1) = TOF(Vdofs,1) + FCbt;
    end
    U = TOK \ TOF;
    SN(:,6) = U(:,1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % for the first macro node of coarse element in the x direction
    % v4
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    TOK = SK;
    TOF = SF;
    
    TOK(s1(ns1,1)*2-1,s1(ns1,1)*2-1) = C*TOK(s1(ns1,1)*2-1,s1(ns1,1)*2-1); 
    TOF(s1(ns1,1)*2-1,1) = 0;
    TOK(s1(ns1,1)*2,s1(ns1,1)*2) = C*TOK(s1(ns1,1)*2,s1(ns1,1)*2);
    TOF(s1(ns1,1)*2,1) = 0;
    
    S3y = linspace(1,0,ns3);
    for i = 1:ns3
        Vdofs = [s3(i,1)*2 s1(i,1)*2];
        TOK(Vdofs,Vdofs) = TOK(Vdofs,Vdofs) + Cbt;
        Udofs = [s3(i,1)*2-1 s1(i,1)*2-1];
        TOK(Udofs,Udofs) = TOK(Udofs,Udofs) + Cbt;
        
        bt0 = S3y(i);
        FCbt = C*[bt0*bt1;bt0*bt2];
        TOF(Vdofs,1) = TOF(Vdofs,1) + FCbt;
    end
    
    S4y = linspace(0,1,ns4);
    for i = 1:(ns4-1)
        Vdofs = [s4(i,1)*2 s2(i,1)*2];
        TOK(Vdofs,Vdofs) = TOK(Vdofs,Vdofs) + Cbt;
        Udofs = [s4(i,1)*2-1 s2(i,1)*2-1];
        TOK(Udofs,Udofs) = TOK(Udofs,Udofs) + Cbt;
        
        bt0 = S4y(i);
        FCbt = C*[bt0*bt1;bt0*bt2];
        TOF(Vdofs,1) = TOF(Vdofs,1) + FCbt;
    end
    U = TOK \ TOF;
    SN(:,8) = U(:,1);
        
    %% Calculate the macroscopic equivalent matrix
    % 
    Q = sparse(2*SNN,1);
    
    % external force is imposed on the right edge
    q = -0.05;
    % X direction
%     Q((2*s2(1)-1),1) = q/(ns2-1)/2;
%     Q((2*s2(ns2)-1),1) = q/(ns2-1)/2;
%     Q((2*s2(2:(ns2-1))-1),1) = q/(ns2-1);
    % Y direction
    Q((2*s2(1)),1) = q/(ns2-1)/2;
    Q((2*s2(ns2)),1) = q/(ns2-1)/2;
    Q((2*s2(2:(ns2-1))),1) = q/(ns2-1);


    % Y direction, external force is imposed on the upper edge
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