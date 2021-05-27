function [ClCdScore, alpha, Cm] = GetClCdScore(Cl,modairfoil_filename)
    
    strtCiter = 50;
    maxCiter = 4*strtCiter;
    Citer = strtCiter;
    full = false;
    bumpCl = 0.02;
    bCount = 0;
    
    while ~full
        ClCdData(Cl,modairfoil_filename,Citer)

         % Extract Cd and alpha corresponding to Cl from setCldata.pol
        FID = fopen('setCldata.pol', 'r');
        clData = textscan(FID,'%f %f %f %f %f %f %f', 'HeaderLines', 12);  %skips headers
        if cellfun(@isempty,clData) %if data is empty, set ClCd and alpha to zero to stop error
            Citer = Citer + strtCiter;
            if Citer >= maxCiter
                bCount = bCount + 1;
                if bCount > 2
                    full = true;
                    alpha = 0;
                    ClCdScore = 0;
                    Cm = -0.01;
                    fprintf('\t\tClCdScore could not generate proper data.\n');
                else
                    Cl = Cl + bumpCl;
                    Citer = strtCiter;
                end
            end
        else
            full = true;
            alpha = clData{1,1}(1);                     %Alpha
            cl = clData{1,2}(1);                        %Coefficient of Lift
            cd = clData{1,3}(1);                        %Coefficient of Drag
            Cm = clData{1,5}(1);                        %Coefficient of Moment
            ClCdScore = -cl/cd;                         %-1*drag efficiency
            fprintf(['\t\tClCdScore = ',num2str(ClCdScore), ', alpha = ', num2str(alpha), ', Cm = ', num2str(Cm),...
                    '\n\t\t\tGenerated data using Cl = ', num2str(Cl), '.\n']);
        end
        fclose('all');
    end
end
    
