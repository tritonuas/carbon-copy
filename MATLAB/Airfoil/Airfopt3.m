%% Major changes from Airfopt2
%Continuous airfoils using Xfoil TSET and HIGH function to set 
%   max camber,max thickness, location of both;
%Functional second-derivative based optimization;
%Stall properties considered;


%% Input Formatting
% Cl = double;
% GuessAirfoil = 'ABCC' using nacaABCC numbering system;
function [best_airfoil, ClCd] = Airfopt3(Cl, GuessAirfoil)
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
    
%% Preallocating arrays?
    %Is this necessary? MATLAB recommends

while converged == false && iter <= maxiter
    %% Record A, B, C, D
    As(iter) = A;
    Bs(iter) = B;
    Cs(iter) = C;
    Ds(iter) = D;
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
%                     
        Airfoils{1,1}(1) = A;      Airfoils{1,1}(2) = B;       Airfoils{1,1}(3) = C;       Airfoils{1,1}(4) = D;      % a
        Airfoils{2,1}(1) = A+dX;   Airfoils{2,1}(2) = B;       Airfoils{2,1}(3) = C;       Airfoils{2,1}(4) = D;      % b
        Airfoils{3,1}(1) = A;      Airfoils{3,1}(2) = B+dX;    Airfoils{3,1}(3) = C;       Airfoils{3,1}(4) = D;      % c
        Airfoils{4,1}(1) = A;      Airfoils{4,1}(2) = B;       Airfoils{4,1}(3) = C+dX;    Airfoils{4,1}(4) = D;      % d
        Airfoils{5,1}(1) = A;      Airfoils{5,1}(2) = B;       Airfoils{5,1}(3) = C;       Airfoils{5,1}(4) = D+dX;   % e
        Airfoils{6,1}(1) = A+2*dX; Airfoils{6,1}(2) = B;       Airfoils{6,1}(3) = C;       Airfoils{6,1}(4) = D;      % f
        Airfoils{7,1}(1) = A+dX;   Airfoils{7,1}(2) = B+dX;    Airfoils{7,1}(3) = C;       Airfoils{7,1}(4) = D;      % g
        Airfoils{8,1}(1) = A+dX;   Airfoils{8,1}(2) = B;       Airfoils{8,1}(3) = C+dX;    Airfoils{8,1}(4) = D;      % h
        Airfoils{9,1}(1) = A+dX;   Airfoils{9,1}(2) = B;       Airfoils{9,1}(3) = C;       Airfoils{9,1}(4) = D+dX;   % i
        Airfoils{10,1}(1) = A;     Airfoils{10,1}(2) = B+2*dX; Airfoils{10,1}(3) = C;      Airfoils{10,1}(4) = D;     % j
        Airfoils{11,1}(1) = A;     Airfoils{11,1}(2) = B+dX;   Airfoils{11,1}(3) = C+dX;   Airfoils{11,1}(4) = D;     % k
        Airfoils{12,1}(1) = A;     Airfoils{12,1}(2) = B+dX;   Airfoils{12,1}(3) = C;      Airfoils{12,1}(4) = D+dX;  % l
        Airfoils{13,1}(1) = A;     Airfoils{13,1}(2) = B;      Airfoils{13,1}(3) = C+2*dX; Airfoils{13,1}(4) = D;     % m
        Airfoils{14,1}(1) = A;     Airfoils{14,1}(2) = B;      Airfoils{14,1}(3) = C+dX;   Airfoils{14,1}(4) = D+dX;  % n
        Airfoils{15,1}(1) = A;     Airfoils{15,1}(2) = B;      Airfoils{15,1}(3) = C;      Airfoils{15,1}(4) = D+2*dX;% o
        
            
            for i = 1:15 %For each airfoil: find score
                %% Get naca airfoil approximating current airfoil (rounding etc.)
                roundA = num2str(round(100*Airfoils{i,1}(1)));
                roundB = num2str(round(10*Airfoils{i,1}(2)));
                roundC = num2str(round(100*Airfoils{i,1}(3)));
                if length(roundC) == 1
                    roundC = ['0' roundC];
                end
                naca = [roundA roundB roundC];
            end
                %% Modify naca airfoil using TSET and HIGH
            fid = fopen('xfoil_input.txt','w');   % create the inputs file 
            fprintf(fid,['naca ', naca, '\n']); %loads rounded airfoil
            fprintf(fid,'gdes\n');   % opens gdes routine
            fprintf(fid,['tset\n' num2str(C) '\n' num2str(A) '\n']); %sets exact max camber and thickness
            fprintf(fid,['high\n' num2str(D) '\n' num2str(B) '\n']); %sets exact location of max c and t
            fprintf(fid,['name\n' 'modairfoil\n']); %renames buffer airfoil
            fprintf(fid,'x\n\n'); %sets modified airfoil as current airfoil
            
            %% Find Cl/Cd Score
            fprintf(fid,'pane\n');   % makes the airfoils nice
            fprintf(fid,'oper\n');fprintf(fid,'visc\n');   % makes visc analysis w/ Re=4E5
            fprintf(fid,'4e5\n');fprintf(fid,'pacc\n');
            fprintf(fid,['setCldata.pol\n\n']);%specifies file to store data
            fprintf(fid,'iter\n');fprintf(fid,'50\n');
            fprintf(fid,['cl ', num2str(Cl),'\n']);% run xfoil with given Cl value
            fprintf(fid,'pacc\n');

            cmd = 'xfoil.exe < xfoil_input.txt';   % running on xfoil
            [~,~] = system(cmd);

            fclose('all');
            delete('xfoil_input.txt');
            
           %% Extract Cd and alpha corresponding to Cl from setCldata.pol
            FID = fopen('setCldata.pol', 'r');
                   D = textscan(fID,'%f %f %f %f %f %f %f', 'HeaderLines', 12);  %skips headers
                    if cellfun(@isempty,D) %if data is empty, set ClCd to zero to stop error
                        alpha = 0;
                        ClCd = 0;
                    else
                        alpha = D{1,1};                     %Alpha
                        cl = D{1,2};                        %Coefficient of Lift
                        cd = D{1,3};                        %Coefficient of Drag
                        ClCd = cl/cd;                       %drag efficiency
                    end
                    fclose('all');
                %Store ClCd in corresponding entry in ClCd  matrix
                %Find stall score   (I:alpha, coordinate filename;
                %                    O:Stall score)
                    %Find max Cl/Cd vs alpha starting at input alpha
                        %
                    %At max Cl/Cd vs alpha, calculate second derivative
                        %
                    %Assign stall score based on second derivative to airfoil
                        % lower |second derivative| --> better: either penalized
                        % less or rewarded more
                        % Weight depends on how strongly optimixation should
                        % value stall properties
                    %Put in corresponding entry in stall score matrix
                % Record Cl/Cd at ABCD: ClCd(1,1)
        
        %Add Cl/Cd score and stall score
        %Find first partial derivatives (fpd) at all necessary points
            %Columns 2:5 - Column 1 of score matrix
            %(Will be 4x5 matrix)
        %Find second partial derivatives (spd) at all necessary points
            %Rows 2:5 - Row 1 of fpd matrix
            %(Will be 4x4 matrix)
    %Use Hessian and Gradient to find Step size
        %Gradient = fpd(1,:)
        %Hessian = spd(:)
        %step = H\G
    %Implement step
        %if updating a variable (v) would break a boundary condition:
            %v = boundary broken    
        %else
            %A = A + step(1)
            %B = B + step(2)
            %C = C + step(3)
            %D = D + step(4)
end
%best airfoil is final recorded airfoil
%Best Cl/Cd corresponds to best airfoil

end





















