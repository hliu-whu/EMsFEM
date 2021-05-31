function OutPut(Nodes,Elems,U,V,filename)

    NN = size(Nodes,1);
    NE = size(Elems,1);
    
    filename = strcat('../res/',filename);
    if(exist(filename,'file')~=0)
        delete(filename);
    end
    
    fid = fopen(filename,'w');

    fprintf(fid,'TITLE = "Results"\n');
    fprintf(fid,'VARIABLES = "X", "Y", "U", "V"\n');
    fprintf(fid,'ZONE N= %8d, E= %8d, F=FEPOINT, ET=QUADRILATERAL\n',NN,NE);
    
    fprintf(fid,'%25.15E %25.15E %25.15E %25.15E\n',[Nodes(:,1) Nodes(:,2) U V]');
    fprintf(fid,'%10d %10d %10d %10d\n',[Elems(:,1) Elems(:,2) Elems(:,3) Elems(:,4)]');
    
    fclose(fid);



end