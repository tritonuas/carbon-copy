function ClCd = findclcd2(Cl, Airfoil, airfoil_dir_name)
%   airfoil_dir_name = 'airfoil_database\\specific_cl\\';

    fID = fopen(strcat(airfoil_dir_name, 'cl', num2str(Cl),'naca', Airfoil, '.pol'), 'r');
    if fID == -1 % If file does not exist, generate file
        fID = genXfoil2(Cl,Airfoil,airfoil_dir_name);
    end
    D = textscan(fID,'%f %f %f %f %f %f %f', 'HeaderLines', 12);  %skips headers
    if cellfun(@isempty,D) %if data is empty, Xfoil didn't generate it
        fprintf('\n Problem Airfoil: %s. Xfoil did not converge.\n',Airfoil);
    end
    fclose(fID);

%   alpha = D{1,1};                     %Alpha
    cl = D{1,2};                        %Coefficient of Lift
    cd = D{1,3};                        %Coefficient of Drag
    ClCd = cl/cd;                       %drag efficiency

    fclose('all');
end