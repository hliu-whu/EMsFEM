function [D] = Get_D(E0,emu0)
% calculate the elastic matrix

    D = zeros(3,3);
    
    % plane strain
%     E0 = E0/(1-emu0^2);
%     emu0 = emu0/(1-emu0);
%     D0 = E0/(1-emu0*emu0);

    % plane stress
    D0 = E0/(1-emu0*emu0);

    
    D(1,1) = 1;
    D(1,2) = emu0;
    D(2,1) = emu0;
    D(2,2) = 1;
    D(3,3) = (1-emu0)/2;
    
    D = D0*D;
    
end