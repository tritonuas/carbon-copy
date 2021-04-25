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
    A = str2double(GuessAirfoil(1))/100;% max camber: chord ratio
    B = str2double(GuessAirfoil(2))/10;% location of max camber : chord ratio
    C = str2double(GuessAirfoil(3:4))/100;% max thickness : Chord ratio
    D = 0.3; % loc of max thickness : chord ratio, always 0.3 for naca numbering
    
    dX = 0.01; %in derivative calculations, (x2-x1)
    Amin = 0;
    Amax = 0.1 - 3*dX; %allows calculation of second derivative
    Bmin = 0;          
    Bmax = 1 - 3*dX;
    Cmin = 0;
    Cmax = 1 - 3*dX;
    Dmin = 0;
    Dmax = 1 - 3*dX;
    
    iter = 1; %iteration counter
    maxiter = 100; %shouldn't be necessary
    
    converged = false; %to enter while loop
    con = 0.001;% converged if max(step) below this value
    modairfoil_filename = 'C:\\TUAS\\carbon-copy\\MATLAB\\Airfoil\\modairfoil.dat';
    
    %% Preallocating arrays
    As = NaN*zeros(maxiter,1);
    Bs = NaN*zeros(maxiter,1);
    Cs = NaN*zeros(maxiter,1);
    Ds = NaN*zeros(maxiter,1);
    ClCds = NaN*zeros(maxiter,1);
    while ~converged && iter <= maxiter
        %% Record A, B, C, D
        As(iter,1) = A;
        Bs(iter,1) = B;
        Cs(iter,1) = C;
        Ds(iter,1) = D;
        
        %% Find Hessian and Gradient
            %Create matrix of all airfoils necessary:
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
        
        Airfoils = [A B C D;...         %a
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
        
        % Preallocating Arrays
        ClCdScore = NaN*zeros(15,1);
        StallScore = NaN*zeros(15,1);
        F = NaN*zeros(15,1);
        fpds = NaN*zeros(5,4);
        Gradient = NaN*zeros(4,1);
        Hessian = NaN*zeros(4,4);
        
        for i = 1:15 %For each airfoil: find score
            modify_airfoil(Airfoils(i,:),modairfoil_filename); % Create coordinates of airfoil

            [ClCdScore(i), alpha] = GetClCdScore(Cl,modairfoil_filename); % Find ClCd Score and starting alpha

            [StallScore(i), alphas, Aclcds] = GetStallScore(alpha,modairfoil_filename); % Find stall score 

            % Plot alpha optimization to check

            StallWeight = 1;
            F(i) = ClCdScore(i) + StallWeight*StallScore(i); %Add Cl/Cd score and stall score
        end
        ClCds(iter) = ClCdScore(1); %ClCd at ABCD
                
%                     a b c d e   
%                     b f g h i
%                     c g j k l
%                     d h k m n
%                     e i l n o
                
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
        Func(4,5) = F(14); Func(5,4) = F(14);   %n
        Func(5,5) = F(15);                      %o

                
        for c = 2:5
            fpds(:,c-1) = (Func(:,c) - Func(:,1))/dX;
                % first row is first partial derivatives at ABCD
                % second row is first partial derivatives at (A+1)BCD
                % third row is first partial derivatives at A(B+1)CD
                % fourth row is first partial derivatives at AB(C+1)D
                % fifth row is first partial derivatives at ABC(D+1)
        end
        for i = 1:4
            % Gradient is column matrix of all first partial derivatives at
            % Airfoil in (1,1)
            Gradient(i,1) = fpds(1,i);
        end
        for r = 2:5
        % Subtract first row from every other row
        % (y2-y1)/(x2-x1), but (x2-x1) always = 1
        Hessian(r-1,:) = (fpds(r,:) - fpds(1,:))/dX;
        end

        step = Hessian\Gradient; %Positive for now because using positive clcd values


        %% Implement step
        if step < con
            converged = true;
        else
            if A + step(1) > Amax
                A = Amax;
            elseif A + step(1) < Amin
                A = Amin;
            else
                A = A + step(1);
            end

            if B + step(2) > Bmax
                B = Bmax;
            elseif B + step(2) < Bmin
                B = Bmin;
            else
                B = B + step(2);
            end

            if C + step(3) > Cmax
                C = Cmax;
            elseif C + step(3) < Cmin
                C = Cmin;
            else
                C = C + step(3);
            end

            if D + step(4) > Dmax
                D = Dmax;
            elseif D + step(4) < Dmin
                D = Dmin;
            else
                D = D + step(4);
            end
        end
    end
best_airfoil_specs = [As(end) Bs(end) Cs(end) Ds(end)];
ClCd = ClCds(end);


end

