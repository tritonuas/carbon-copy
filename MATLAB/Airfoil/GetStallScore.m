function [StallScore, alphas, aclcds] = GetStallScore(alpha,modairfoil_filename)
 
    Aconverged = false;
    Aiter = 0;
    inc = 0.5; % delta alpha for Xfoil
    Acon = 0.02;% Converged once Astep below this 
    Amaxiter = 50;
    maxStep = 5;
    strtPanel = 160;
    maxPanel = 400;
    Piter = 40;

    while ~Aconverged && Aiter < Amaxiter% Find max Cl/Cd vs alpha starting at input alpha
        Aiter = Aiter + 1;
        alphas(Aiter) = alpha;
        
        panel = strtPanel;
        full = false;
        
        while full == false
            AlphaData2(alpha, modairfoil_filename, inc, panel)   %runs Xfoil for necessary alphas, stores data in 
                                                                %alphadata.pol
            FID = fopen('alphadata.pol', 'r');
            aData = textscan(FID,'%f %f %f %f %f %f %f', 'HeaderLines', 12);  %skips headers
                acl = aData{:,2};                        %Coefficient of Lift
                acd = aData{:,3};                        %Coefficient of Drag
                aclcd = acl./acd;                           %drag efficiency
            fclose('all');
            
            disp(aclcd);
            
            if length(aclcd) == 3
                full = true;
            else
                panel = panel + Piter; %increase panel number until converges
                if panel >= maxPanel
                    alpha = alpha + inc/2; % Shifts slightly and retries
                    panel = strtPanel;
%                     fprintf(['Xfoil is having trouble converging for the current airfoil.\n',...
%                             'Please use Xfoil to generate Cl and Cd values for the airfoil at:\n',...
%                             '\talpha = ', num2str(alpha), ':', num2str(inc), ':', num2str(alpha+2*inc) '\n',...
%                             'and store them in alphadata.pol.\n\n'...
%                             'Tricks to get Xfoil to converge:\n',...
%                             '\t1) Try messing with paneling size using ppar command.\n',...
%                             '\t2) Xfoil converges better with an open tail. ',...
%                                 'Delete some points from the start and end of ',...
%                                 modairfoil_filename 'and rerun Xfoil.\n\n']);
%                             input('Press any key when correct data generated.    ')
                end
            end
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
        else
            alpha = alpha + aStep;
        end
    end
    StallScore = -1/ddA(aclcds == max(aclcds)); %lower ddA is better, all should be negative
    
%     figure
%     plot(1:Aiter,alphas)
%     xlabel('Iteration'); ylabel('Alpha');
%     
%     figure
%     plot(1:Aiter,aclcds)
%     xlabel('Iteration'); ylabel('Cl/Cd');

%     figure
%     plot(1:Aiter,ddA)
%     xlabel('Iteration'); ylabel('Derivative');
    
end