function FID = genXfoil2(Cl, Airfoil, airfoil_dir_name)
%     airfoil_dir_name = airfoil_database\\specific_cl\\;

    fid = fopen('xfoil_input.txt','w');   % create the inputs file 
    fprintf(fid,['naca ', Airfoil, '\n']);
    fprintf(fid,'pane\n');   % makes the airfoils nice
    fprintf(fid,'oper\n');fprintf(fid,'visc\n');   % makes visc analysis w/ Re=4E5
    fprintf(fid,'4e5\n');fprintf(fid,'pacc\n');
    fprintf(fid,[airfoil_dir_name,'cl', num2str(Cl), 'naca', Airfoil, '.pol\n\n']);
    fprintf(fid,'iter\n');fprintf(fid,'50\n');
    fprintf(fid,['cl ', num2str(Cl),'\n']);% run xfoil with given Cl value
    fprintf(fid,'pacc\n');

    cmd = 'xfoil.exe < xfoil_input.txt';   % running on xfoil
    [~,~] = system(cmd);

    fclose('all');
    delete('xfoil_input.txt');
    FID = fopen([airfoil_dir_name,'cl', num2str(Cl), 'naca', Airfoil, '.pol'], 'r');
end