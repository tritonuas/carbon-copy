function AlphaData2(alpha,modairfoil_filename, inc, panel)
    % Generate and save data for set of alpha
    
    if exist('alphadata.pol','file') == 2
        delete('alphadata.pol');%so Xfoil can write new data there
    end
    fid = fopen('xfoil_input.txt','w');
    fprintf(fid, ['load\n', modairfoil_filename,'\n']);% loads modified airfoil coordinates
    fprintf(fid, ['ppar\n','n\n', num2str(panel), '\n\n\n']);
    fprintf(fid,'oper\n');
    fprintf(fid,['visc\n', '4e5\n']);  % makes visc analysis w/ Re=4E5
    fprintf(fid,'pacc\n');
    fprintf(fid,'alphadata.pol\n\n');%specifies file to store data
    fprintf(fid,['iter\n', '100\n']); % set iteration number to 50
    fprintf(fid,'aseq\n');% run xfoil with alpha sequence
    fprintf(fid,[num2str(alpha), '\n', num2str(alpha+2*inc), '\n', num2str(inc), '\n']);
                          %aseq format:first alpha, last alpha, increment size
    fprintf(fid,'pacc\n');
    cmd = 'xfoil.exe < xfoil_input.txt';   % running on xfoil
    [~,~] = system(cmd);

    fclose('all');
    delete('xfoil_input.txt');
end