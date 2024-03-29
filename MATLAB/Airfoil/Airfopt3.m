%% Major changes from Airfopt2
%Continuous airfoils using Xfoil TSET and HIGH function to set 
%   max camber,max thickness, location of both;
%Functional second-derivative based optimization;
%Stall properties considered;

%% Input Formatting
% Cl = double;
% GuessAirfoil = 'ABCC' using nacaABCC numbering system;

function [best_airfoil_specs, ClCd] = Airfopt3(Cl, GuessAirfoil)
%% Set starting values
    dX = 0.01; %in derivative calculations, (x2-x1)
    minA = 0; 
    maxA = 0.1;
    minB = 0.1;          
    maxB = 0.6;
    minC = 0.1;
    maxC = 0.3;
    minD = 0.2;
    maxD = 0.6;
    minCm = -0.1;
    maxCm = 0;
    
    %Constraint Functions
        bowl = @(x,y,z) 1/(0.5*(x-y)) - 1/(0.5*(x-z));%y = min, z = max
    
    A = str2double(GuessAirfoil(1))/100;% max camber: chord ratio
        if A > maxA
            A = maxA;
        elseif A < minA %don't need to check min as long as min == 0
            A = minA;
        end
    B = str2double(GuessAirfoil(2))/10;% location of max camber : chord ratio
        if B > maxB
            B = maxB;
        elseif B < minB %don't need to check min as long as min == 0
            B = minB;
        end
    C = str2double(GuessAirfoil(3:4))/100;% max thickness : Chord ratio
        if C > maxC
            C = maxB;
        elseif C < minC %don't need to check min as long as min == 0
            C = minC;
        end
    D = 0.3; % loc of max thickness : chord ratio, always 0.3 for naca numbering
        if D > maxD %Starting D should always be within range
            D = maxD;
        elseif D < minD %don't need to check min as long as min == 0
            D = minD;
        end     
        
    iter = 0; %iteration counter
    maxiter = 100; %shouldn't be necessary if optimization works properly
    converged = false; %to enter while loop
    con = 0.001;% converged if max(step) below this value
    modairfoil_filename = 'modairfoil.dat';
%   modairfoil_filename = 'C:\\TUAS\\carbon-copy\\MATLAB\\Airfoil\\modairfoil.dat';
%     bestairfoil_filename = 'C:\\TUAS\\carbon-copy\\MATLAB\\Airfoil\\bestairfoil.dat';
    StallWeight = 1;
    
    %% Preallocating arrays
%     As = NaN*zeros(maxiter,1); %not possible if maxiter doesn't exist
%     Bs = NaN*zeros(maxiter,1);
%     Cs = NaN*zeros(maxiter,1);
%     Ds = NaN*zeros(maxiter,1);
%     ClCds = NaN*zeros(maxiter,1);
    
    while ~converged && iter <= maxiter
        iter = iter + 1;
        
        % Record A, B, C, D
        As(iter,1) = A;
        Bs(iter,1) = B;
        Cs(iter,1) = C;
        Ds(iter,1) = D;
        
       bunch = [As(iter,1) Bs(iter,1) Cs(iter,1) Ds(iter,1)];
       fprintf(['Iteration ', num2str(iter),': ' num2str(bunch), '\n']);
            %% Create matrix of all airfoils necessary:
                %Will be cell array of double arrays
                    %  ABCD      (A+dX)BCD       A(B+dX)CD       AB(C+dX)D       ABC(D+dX)
                    % (A+dX)BCD  (A+2dX)BCD     (A+dX)(B+dX)CD  (A+dX)B(C+dX)D  (A+dX)BC(D+dX)
                    %  A(B+dX)CD (A+dX)(B+dX)CD  A(B+2dX)CD      A(B+dX)(C+dX)D  A(B+dX)C(D+dX)
                    %  AB(C+dX)D (A+dX)B(C+dX)D  A(B+dX)(C+dX)D  AB(C+2dX)D      ABC(C+dX)(D+dX) 
                    %  ABC(D+dX) (A+dX)BC(D+dX)  A(B+dX)C(D+dX)  AB(C+dX)(D+dX)  ABC(D+2dX)

%                     a b c d e   
%                     b f g h i
%                     c g j k l
%                     d h k m n
%                     e i l n o
        
        Airfoil_Specs = [A B C D;...         %a
                         A+dX B C D;...      %b
                         A B+dX C D;...      %c
                         A B C+dX D;...      %d
                         A B C D+dX;...      %e
                         A+2*dX B C D;...    %f
                         A+dX B+dX C D;...   %g
                         A+dX B C+dX D;...   %h
                         A+dX B C D+dX;...   %i
                         A B+2*dX C D;...    %j
                         A B+dX C+dX D;...   %k
                         A B+dX C D+dX;...   %l
                         A B C+2*dX D;...    %m
                         A B C+dX D+dX;...   %n
                         A B C D+2*dX];      %o
                                                
        %% Generate ClCd and Stall Propterty Data
        % Preallocating Arrays
        ClCdScore = NaN*zeros(15,1);
        StallScore = NaN*zeros(15,1);
        F = NaN*zeros(15,1);
        fpds = NaN*zeros(5,4);
%         Gradient = NaN*zeros(4,1);
        Hessian = NaN*zeros(4,4);
        
        for i = 1:15 %For each airfoil: find score
%             disp(Airfoil_Specs(i,:)); %TESTING
            fprintf(['\tTesting Airfoil: ', num2str(Airfoil_Specs(i,:)), '\n']);            
            
            tic
            modify_airfoil(Airfoil_Specs(i,:),modairfoil_filename); % Create coordinates of airfoil
            modtime(i,iter) = toc;
%             if i == 1 %This airfoil is better than last iteration, Specs(1,:) is airfoil being analyzed
%                 %delete(bestairfoil_filename);
%                 [success, message]= copyfile(modairfoil_filename,bestairfoil_filename);
% %                 disp(success); disp(message);%TESTING
%             end
            
            tic
            [ClCdScore(i,1), alpha, Cm(i,1)] = GetClCdScore(Cl,modairfoil_filename); % Find ClCd Score and starting alpha
            clcdtime(i,iter) = toc;
            
            tic
% %             [StallScore(i,1),~,~] = GetStallScore(alpha,modairfoil_filename); % Find stall score 
%             [StallScore(i),~,~] = BFGetStallScore(modairfoil_filename);
            StallScore(i,1) = 0;
            stalltime(i,iter) = toc;
            
            % Objective Function
            F(i) = ClCdScore(i,1) + StallWeight*StallScore(i,1)...
                    + bowl(Airfoil_Specs(i,1), minA, maxA) + bowl(Airfoil_Specs(i,2), minB, maxB)...
                    + bowl(Airfoil_Specs(i,3), minC, maxC) + bowl(Airfoil_Specs(i,4), minD, maxD)...
                    + bowl(Cm(i,1), minCm, maxCm);
            
            fprintf(['\t\tObjective Score: ', num2str(F(i)), '\n']);
            fprintf(['\t\tTiming: modairfoil = ', num2str(modtime(i,iter)), 'sec, ClCdScore = ', num2str(clcdtime(i,iter)), 'sec, StallScore = ', num2str(stalltime(i,iter)), 'sec\n']);
        end
        
        ClCds(iter) = ClCdScore(1); %ClCd at ABCD
        Objective(iter) = F(1);    
        
%         Assemble matrix:        
%             a b c d e   
%             b f g h i
%             c g j k l
%             d h k m n
%             e i l n o

        Func(1,1) = F(1);                       %a
        Func(2,1) = F(2); Func(1,2) = F(2);     %b
        Func(3,1) = F(3); Func(1,3) = F(3);     %c
        Func(4,1) = F(4); Func(1,4) = F(4);     %d
        Func(5,1) = F(5); Func(1,5) = F(5);     %e
        Func(2,2) = F(6);                       %f
        Func(3,2) = F(7); Func(2,3) = F(7);     %g
        Func(4,2) = F(8); Func(2,4) = F(8);     %h
        Func(5,2) = F(9); Func(2,5) = F(9);     %i
        Func(3,3) = F(10);                      %j
        Func(4,3) = F(11); Func(3,4) = F(11);   %k
        Func(5,3) = F(12); Func(3,5) = F(12);   %l
        Func(4,4) = F(13);                      %m
        Func(5,4) = F(14); Func(4,5) = F(14);   %n
        Func(5,5) = F(15);                      %o

       %% Derivative and Step Calculations         
        
       fpds(:,1:4) = (Func(:,2:5)-Func(:,1))/dX;
            % first row is first partial derivatives at ABCD
            % second row is first partial derivatives at (A+1)BCD
            % third row is first partial derivatives at A(B+1)CD
            % fourth row is first partial derivatives at AB(C+1)D
            % fifth row is first partial derivatives at ABC(D+1)

        Gradient(:,iter) = fpds(1,:);
            % Gradient is column matrix of all first partial derivatives at
            % Airfoil in (1,1)
            
        Hessian(1:4,:) = (fpds(2:5,:) - fpds(1,:))/dX;
            % Subtract first row from every other row
            % (y2-y1)/(x2-x1)
       
        step = abs(Hessian)\-Gradient(:,iter);
        %% Implement step
        if max(abs(step)) < con
            converged = true;
        else
            A = A + step(1);
%             if A > maxA
%                 A = maxA;
%             elseif A < minA
%                 A = minA;
%             end
            
            B = B + step(2);
%             if B > maxB
%                 B = maxB;
%             elseif B < minB
%                 B = minB;
%             end
            
            C = C + step(3);
%             if C > maxC
%                 C = maxC;
%             elseif C < minC
%                 C = minC;
%             end
            
            D = D + step(4);
%             if D > maxD
%                 D = maxD;
%             elseif D < minD
%                 D = minD;
%             end
        end
    end
best_airfoil_specs = [As(end) Bs(end) Cs(end) Ds(end)];
ClCd = ClCds(end);

figure
plot(1:iter,ClCds);
xlabel('Iteration');
ylabel('Cl/Cd');

figure
plot(1:iter,As);
xlabel('Iteration');
ylabel('Max Camber');

figure
plot(1:iter,Bs);
xlabel('Iteration');
ylabel('Location of Max Camber');

figure
plot(1:iter,Cs);
xlabel('Iteration');
ylabel('Max Thickness');

figure
plot(1:iter,Ds);
xlabel('Iteration');
ylabel('Location of Max Thickness');

disp([As Bs Cs Ds]);

figure
plot(1:iter,Gradient(1,:));
xlabel('Iteration');
ylabel('A Gradient');

figure
plot(1:iter,Gradient(2,:));
xlabel('Iteration');
ylabel('B Gradient');

figure
plot(1:iter,Gradient(3,:));
xlabel('Iteration');
ylabel('C Gradient');

figure
plot(1:iter,Gradient(4,:));
xlabel('Iteration');
ylabel('D Gradient');

figure
plot(1:iter,Objective);
xlabel('Iteration');
ylabel('Objective');

end

