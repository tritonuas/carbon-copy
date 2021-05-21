function [StallScore, maxAlpha, maxClCd] = BFGetStallScore(modairfoil_filename)
    if exist('alphadata.pol','file') == 2
        delete('alphadata.pol');%so Xfoil can write new data there
    end
    
    inc = 0.5; %increment size
    tol = inc/10; %tolerance to make == not act up
    
    fid = fopen('xfoil_input.txt','w');
    fprintf(fid, ['load\n', modairfoil_filename,'\n']);% loads modified airfoil coordinates
    fprintf(fid,'pane\n');   % makes the airfoils nice
    fprintf(fid,'oper\n');
    fprintf(fid,['visc\n', '4e5\n']);  % makes visc analysis w/ Re=4E5
    fprintf(fid,'pacc\n');
    fprintf(fid,'alphadata.pol\n\n');%specifies file to store data
    fprintf(fid,['iter\n', '50\n']); % set iteration number to 50
    fprintf(fid,'aseq\n');% run xfoil with alpha sequence
    fprintf(fid,['0\n20\n',num2str(inc),'\n']);
                          %aseq format:first alpha, last alpha, increment size
    fprintf(fid,'pacc\n');

    cmd = 'xfoil.exe < xfoil_input.txt';   % running on xfoil
    [~,~] = system(cmd);
    
    fclose('all');
    delete('xfoil_input.txt');
    
    fID = fopen('alphadata.pol','r');  %loads the file
    AlphaData = textscan(fID,'%f %f %f %f %f %f %f', 'HeaderLines', 12);  %skips headers
    fclose(fID);
    alpha = AlphaData{:,1};                     %range of alpha
    cl = AlphaData{:,2};                        %Coefficient of Lift
    cd = AlphaData{:,3};                        %Coefficient of Drag
    aclcd = cl./cd;
    
    maxClCd = max(aclcd);
    maxAlpha = alpha(aclcd == maxClCd);
    dA1 = (aclcd(abs(maxAlpha+inc-alpha) < tol)-aclcd(abs(maxAlpha-alpha) < tol))/inc; %(y2-y1)/(x2-x1) first derivative
    dA2 = (aclcd(abs(maxAlpha+2*inc-alpha) < tol)-aclcd(abs(maxAlpha+inc-alpha) < tol))/inc; %(y3-y2)/(x3-x2)
    ddA = (dA2-dA1)/inc;%second derivative
    StallScore = -1/ddA;%change to reflect larger trend

%     figure;
%     hold on
%     plot(alpha, aclcd);
%     xlabel('alpha');
%     ylabel('Cl/Cd');
%     plot(maxAlpha, maxClCd, 'o');
%     hold off
    
end