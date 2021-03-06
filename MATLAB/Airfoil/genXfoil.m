function FID = genXfoil(directory, Airfoil)
    airfoil_d_name = "airfoil_database";
    addpath(airfoil_d_name);

    fid = fopen('xfoil_input.txt','w');   % create the inputs  
    fprintf(fid,['naca ', Airfoil, '\n']);
    fprintf(fid,'pane\n');   % makes the airfoils nice
    fprintf(fid,'oper\n');fprintf(fid,'visc\n');   % makes visc analysis w/ Re=4E5
    fprintf(fid,'4e5\n');fprintf(fid,'pacc\n');
    fprintf(fid,"airfoil_database/naca"+Airfoil+".pol\n\n");
    fprintf(fid,'iter\n');fprintf(fid,'25\n'); % aint nobody got time fo 50
    fprintf(fid,'aseq 0 20 0.5\n');% alpha 0 to 20 deg in 0.5 deg increment
    fprintf(fid,'pacc\n');

    cmd = 'xfoil.exe < xfoil_input.txt';   % running on xfoil
    [~,~] = system(cmd);

    fclose('all');
    delete('xfoil_input.txt');
    FID = fopen(strcat(directory, 'naca', Airfoil, '.pol'), 'r');
end