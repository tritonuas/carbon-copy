function [StallScore, alphas, Naclcds] = GetStallScore(alpha,modairfoil_filename)
 
    Aconverged = false;
    Aiter = 0;
    inc = 0.5; % delta alpha for Xfoil
    Acon = 0.02;% Converged once Astep below this 
    Amaxiter = 50;% Max iterations
    maxStep = 7;% Funky data makes constraint necessary
    strtPanel = 160; %Panels of airfoil in Xfoil
    maxPanel = 400;
    Piter = 40; % Increase panels by this much
    adjust = 0;% Number of times 3 points not generated

    while ~Aconverged && Aiter < Amaxiter% Find max Cl/Cd vs alpha starting at input alpha
        Aiter = Aiter + 1;
        alphas(Aiter) = alpha;
        
        panel = strtPanel;
        full = false;
        
        while ~full
            AlphaData2(alpha, modairfoil_filename, inc, panel)   %runs Xfoil for necessary alphas, stores data in 
                                                                %alphadata.pol
            FID = fopen('alphadata.pol', 'r');
            aData = textscan(FID,'%f %f %f %f %f %f %f', 'HeaderLines', 12);  %skips headers
                acl = aData{:,2};                           %Coefficient of Lift
                acd = aData{:,3};                           %Coefficient of Drag
                Naclcd = -acl./acd;                         % -1*drag efficiency
            fclose('all');
            
%             disp(aclcd);
            
            if length(Naclcd) == 3
                full = true;
                Adjust(Aiter) = adjust;
                Panel(Aiter) = panel;
            else
                panel = panel + Piter; %increase panel number until converges
                if panel >= maxPanel % If not converging
                    adjust = adjust + 1;
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

        Naclcds(Aiter) = Naclcd(1); %for plotting

        % Derivatives and step calculations
        dA1 = (Naclcd(2)-Naclcd(1))/inc; %(y2-y1)/(x2-x1) first derivative
        dA2 = (Naclcd(3)-Naclcd(2))/inc; %(y3-y2)/(x3-x2)
        ddA(Aiter) = (dA2-dA1)/inc;%second derivative
        aStep = -dA1/abs(ddA(Aiter));
        
        if aStep > maxStep %funky data makes constraints necessary
            aStep = maxStep;
        elseif aStep < -maxStep
            aStep = -maxStep;
        end
        
        if abs(aStep) < Acon
            Aconverged = true;
        elseif Aiter > 1 && any(Naclcds(Aiter)==Naclcds([1:Aiter-1,Aiter+1:end]))%If bouncing between multiple alphas
            Aconverged = true;
        else
            alpha = alpha + aStep;
        end
    end
    
    avgPanel = sum(Panel)/length(Panel);
    avgAdjust = sum(Adjust)/length(Adjust);
    if Aconverged
        fprintf(['\t\tGetStallScore converged in ', num2str(Aiter), ' iterations.\n'...
            '\t\t\tAverage # panels: ' num2str(avgPanel), '   Average alpha adjustments: ', num2str(avgAdjust), '\n']);
    else
        fprintf(['\tGetStallScore did not converge in ', num2str(Amaxiter), ' iterations.\n'...
            '\t\t\tAverage # panels: ' num2str(avgPanel), '   Average alpha adjustments: ', num2str(avgAdjust), '\n']);
    end
%     SS = ddA(Naclcds == min(Naclcds));
%     SS = SS(1);
    SS = sum(ddA)/length(ddA);
    StallScore = 1/SS; %lower ddA is better, all should be negative
        
%     figure
%     plot(1:Aiter,alphas)
%     xlabel('Iteration'); ylabel('Alpha');
%     
%     figure
%     plot(1:Aiter,Naclcds)
%     xlabel('Iteration'); ylabel('Cl/Cd');

%     figure
%     plot(1:Aiter,ddA)
%     xlabel('Iteration'); ylabel('Derivative');
    
end