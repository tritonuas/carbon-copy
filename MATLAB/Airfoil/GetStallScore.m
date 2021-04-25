function [StallScore, alphas, Aclcds] = GetStallScore(alpha,modairfoil_filename)
 
    Aconverged = false;
    Aiter = 1;
    inc = 0.1; % delta alpha
    Acon = 0.01;% Converged once Astep below this value

    while ~Aconverged % Find max Cl/Cd vs alpha starting at input alpha
        alphas(Aiter) = alpha;

        AlphaData(alpha,modairfoil_filename, inc) %runs Xfoil for necessary alphas, stores data in 
                                                  %alphadata.pol

        % Extract data from alphadata.pol
        FID = fopen('alphadata.pol', 'r');
        aData = textscan(FID,'%f %f %f %f %f %f %f', 'HeaderLines', 12);  %skips headers
        %if cellfun(@isempty,aData) %if data is empty, set ClCd to zero to stop error
        %    ClCd = 0;
        %else
            cl = aData{:,2};                        %Coefficient of Lift
            cd = aData{:,3};                        %Coefficient of Drag
            ClCd = cl/cd;                           %drag efficiency
        %end

        fclose('all');

        Aclcds(Aiter) = ClCd(1); %for plotting

        % Derivatives and step calculations
        dA1 = (ClCd(2)-ClCd(1))/inc; %(y2-y1)/(x2-x1)
        dA2 = (ClCd(3)-ClCd(2))/inc; %(y3-y2)/(x3-x2)
        ddA = (dA2-dA1)/inc;
        aStep = dA1/ddA;
        
        if aStep < Acon
            Aconverged = true;
            StallScore = ddA;%Might need to recalculate second derivative for more accurate score
        else
            alpha = alpha + aStep;
            Aiter = Aiter + 1;
        end
    end
end