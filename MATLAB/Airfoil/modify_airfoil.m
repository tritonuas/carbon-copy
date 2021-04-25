function modify_airfoil(Airfoil_Specs, modairfoil_filename)
    % Get naca airfoil approximating current airfoil (rounding etc.)
        roundA = num2str(round(100*Airfoil_Specs(1)));
        roundB = num2str(round(10*Airfoil_Specs(2)));
        roundC = num2str(round(100*Airfoil_Specs(3)));
        if length(roundC) == 1
            roundC = ['0' roundC];
        end
        naca = [roundA roundB roundC];

    % Modify naca airfoil using TSET and HIGH
        fid = fopen('xfoil_input.txt','w');   % create the inputs file 
        fprintf(fid,['naca ', naca, '\n']); %loads rounded airfoil
        fprintf(fid,'gdes\n');   % opens gdes routine
        fprintf(fid,['tset\n' num2str(Airfoil_Specs(3)) '\n' num2str(Airfoil_Specs(1)) '\n']); %sets exact max camber and thickness
        fprintf(fid,['high\n' num2str(Airfoil_Specs(4)) '\n' num2str(Airfoil_Specs(2)) '\n']); %sets exact location of max c and t
        fprintf(fid,['name\n' 'modairfoil\n']); %renames buffer airfoil
        fprintf(fid,'x\n\n'); %sets modified airfoil as current airfoil
        fprintf(fid,['save\n', modairfoil_filename, '\n']);

        cmd = 'xfoil.exe < xfoil_input.txt';   % running on xfoil
        [~,~] = system(cmd);

        fclose('all');
        delete('xfoil_input.txt');
end