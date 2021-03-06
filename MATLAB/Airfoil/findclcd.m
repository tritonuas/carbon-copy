function nClCd = findclcd(CL, Airfoil, fdirectory)
    fID = fopen(strcat(fdirectory, 'naca', Airfoil, '.pol'), 'r');
    if fID == -1 % If file does not exist, generate file
        fID = genXfoil(fdirectory, Airfoil);
    end
    D = textscan(fID,'%f %f %f %f %f %f %f', 'HeaderLines', 12);  %skips headers
    fclose(fID);

    %         alpha = D{:,1};                     %Alpha
    cl = D{:,2};                        %Coefficient of Lift
    cd = D{:,3};                        %Coefficient of Drag
    clcd = cl./cd;                      %drag efficiency

    error = abs(cl-CL);                 % What if there are two places where cl == Cl
    Close = cl(error == min(error));
     
%     disp(Airfoil);% TEMPORARY FOR TESTING
%     disp(Close);% TEMPORARY FOR TESTING
     
    Close = Close(1);

    fclose('all');

    nClCd = -clcd(cl == Close); % -Cl/Cd at Cl Approximation
    
    
% delete([fdirectory, 'naca', Airfoil, '.pol']);% TEMPORARY FOR TESTING

end