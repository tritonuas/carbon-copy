function  ClCd = findclcd(Cl, Airfoil, airfoil_dir_name)
    fID = fopen(strcat(airfoil_dir_name, 'naca', Airfoil, '.pol'), 'r');
    if fID == -1 % If file does not exist, generate file
        fID = genXfoil(Airfoil,airfoil_dir_name);
    end
    D = textscan(fID,'%f %f %f %f %f %f %f', 'HeaderLines', 12);  %skips headers
    fclose(fID);

%   alpha = D{:,1};                     %Alpha
    cl = D{:,2};                        %Coefficient of Lift
    cd = D{:,3};                        %Coefficient of Drag
    clcd = cl./cd;                      %drag efficiency

    error = abs(cl-Cl);                
    Close = cl(error == min(error));    
    Close = Close(1);
    
    ClCd = clcd(cl == Close);           % -Cl/Cd at Cl Approximation
    ClCd = ClCd(1);

    fclose('all');
end