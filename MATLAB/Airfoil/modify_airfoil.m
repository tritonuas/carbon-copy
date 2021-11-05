% xfoil doesn't write over old modairfoil.dat with new coordinates
% Better to delete xfoilinput.txt or just clear contents? Same for modairfoil.dat
function modify_airfoil(Airfoil_Specs, modairfoil_filename)
    %% Get naca airfoil approximating current airfoil (rounding etc.)
        roundA = num2str(round(100*Airfoil_Specs(1)));
            if length(roundA) > 1
                roundA = '9';
            end
        roundB = num2str(round(10*Airfoil_Specs(2)));
            if length(roundB) > 1
                roundB = '9';
            end
        roundC = num2str(round(100*Airfoil_Specs(3)));
            if length(roundC) == 1
                roundC = ['0' roundC];
            elseif length(roundC) > 2
                roundC = '99';
            end
        naca = [roundA roundB roundC];

    %% Modify naca airfoil using TSET and HIGH
        if exist(modairfoil_filename,'file') == 2
            delete(modairfoil_filename); %so Xfoil can write new data there
        end
        fid = fopen('xfoil_input.txt','w');   % create the inputs file 
        fprintf(fid,['naca ', naca, '\n']); %loads rounded airfoil
        fprintf(fid,'gdes\n');   % opens gdes routine
        fprintf(fid,['tset\n' num2str(Airfoil_Specs(3)) '\n' num2str(Airfoil_Specs(1)) '\n']); %sets exact max camber and thickness
        fprintf(fid,['high\n' num2str(Airfoil_Specs(4)) '\n' num2str(Airfoil_Specs(2)) '\n']); %sets exact location of max c and t
        fprintf(fid,['name\n' 'modairfoil\n']); %renames buffer airfoil
        fprintf(fid,'x\n\n'); %sets modified airfoil as current airfoil
        fprintf(fid,['save\n', modairfoil_filename, '\n']);

        % WHEN RUNNING ON WINDOWS
        cmd = 'xfoil.exe < xfoil_input.txt';   % running on xfoil
        % WHEN RUNNING ON LINUX
        %cmd = 'xfoil < xfoil_input.txt';   % running on xfoil
        [~,~] = system(cmd);
        
        fclose('all');
        delete('xfoil_input.txt');
end