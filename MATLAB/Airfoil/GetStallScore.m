function [StallScore, alphas, aclcds] = GetStallScore(alpha,modairfoil_filename)
 
    Aconverged = false;
    Aiter = 1;
    inc = 0.5; % delta alpha for Xfoil
    Acon = 0.01;% Converged once Astep below this 
    Amaxiter = 50;
    maxStep = 5;

    while ~Aconverged && Aiter < Amaxiter% Find max Cl/Cd vs alpha starting at input alpha
        
        alphas(Aiter) = alpha;
        
        AlphaData(alpha,modairfoil_filename, inc) %runs Xfoil for necessary alphas, stores data in 
                                                  %alphadata.pol

        % Extract data from alphadata.pol
        FID = fopen('alphadata.pol', 'r');
        aData = textscan(FID,'%f %f %f %f %f %f %f', 'HeaderLines', 12);  %skips headers
        %if cellfun(@isempty,aData) %if data is empty, set ClCd to zero to stop error
        %    aclcd = 0;
        %else
            acl = aData{:,2};                        %Coefficient of Lift
            acd = aData{:,3};                        %Coefficient of Drag
            aclcd = acl./acd;                           %drag efficiency
        %end
disp(aclcd);
if length(aclcd) < 3
    pause
end
        fclose('all');

        aclcds(Aiter) = aclcd(1); %for plotting

        % Derivatives and step calculations
        dA1 = (aclcd(2)-aclcd(1))/inc; %(y2-y1)/(x2-x1) first derivative
        dA2 = (aclcd(3)-aclcd(2))/inc; %(y3-y2)/(x3-x2)
        ddA(Aiter) = (dA2-dA1)/inc;%second derivative
        aStep = dA1/abs(ddA(Aiter));
        if aStep > maxStep
            aStep = maxStep;
        elseif aStep < -maxStep
            aStep = -maxStep;
        end
        
        if abs(aStep) < Acon
            Aconverged = true;
            StallScore = -1/ddA(Aiter);%lower ddA is better, all should be negative
        elseif Aiter == Amaxiter
            StallScore = -1/ddA(clcds == max(clcds));
        else
            alpha = alpha + aStep;
            Aiter = Aiter + 1;
        end
    end
    
    figure
    plot(1:Aiter,alphas)
    xlabel('Iteration'); ylabel('Alpha');
    
    figure
    plot(1:Aiter,aclcds)
    xlabel('Iteration'); ylabel('Cl/Cd');
    
end