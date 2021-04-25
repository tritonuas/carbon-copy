function [ClCdScore, alpha] = GetClCdScore(Cl,modairfoil_filename)
    % Generate and save data for Cl
        fid = fopen('xfoil_input.txt','w');
        fprintf(fid,['load\n', modairfoil_filename, '\n']);
        fprintf(fid,'pane\n');   % makes the airfoils nice
        fprintf(fid,'oper\n');
        fprintf(fid,['visc\n', '4e5\n']);  % makes visc analysis w/ Re=4E5
        fprintf(fid,'pacc\n');
        fprintf(fid,'setCldata.pol\n\n');%specifies file to store data
        fprintf(fid,['iter\n', '50\n']); % set iteration number to 50
        fprintf(fid,['cl ', num2str(Cl),'\n']);% run xfoil with given Cl value
        fprintf(fid,'pacc\n');

        cmd = 'xfoil.exe < xfoil_input.txt';   % running on xfoil
        [~,~] = system(cmd);

        fclose('all');
        delete('xfoil_input.txt');

    % Extract Cd and alpha corresponding to Cl from setCldata.pol
        FID = fopen('setCldata.pol', 'r');
        clData = textscan(FID,'%f %f %f %f %f %f %f', 'HeaderLines', 12);  %skips headers
        if cellfun(@isempty,clData) %if data is empty, set ClCd and alpha to zero to stop error
            alpha = 0;
            ClCdScore = 0;
        else
            alpha = clData{1,1};                     %Alpha
            cl = clData{1,2};                        %Coefficient of Lift
            cd = clData{1,3};                        %Coefficient of Drag
            ClCdScore = cl/cd;                       %drag efficiency
        end
        fclose('all');
end