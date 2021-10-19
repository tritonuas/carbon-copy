function ClCdData(Cl,modairfoil_filename,Citer)
    % Generate and save data for Cl
        if exist('setCldata.pol','file') == 2
            delete('setCldata.pol'); %so Xfoil can write new data there
        end
        fid = fopen('xfoil_input.txt','w');
        fprintf(fid,['load\n', modairfoil_filename, '\n']);
        fprintf(fid,'pane\n');   % makes the airfoils nice
        fprintf(fid,'oper\n');
        fprintf(fid,['visc\n', '4e5\n']);  % makes visc analysis w/ Re=4E5
        fprintf(fid,'pacc\n');
        fprintf(fid,'setCldata.pol\n\n');%specifies file to store data
        fprintf(fid,['iter\n', num2str(Citer),'\n']); % set iteration number to 50
        fprintf(fid,['cl ', num2str(Cl),'\n']);% run xfoil with given Cl value
        fprintf(fid,'pacc\n');

        % WHEN RUNNING ON WINDOWS
        % cmd = 'xfoil.exe < xfoil_input.txt';   % running on xfoil
        % WHEN RUNNING ON LINUX
        cmd = 'xfoil < xfoil_input.txt';   % running on xfoil
        [~,~] = system(cmd);

        fclose('all');
        delete('xfoil_input.txt');
end