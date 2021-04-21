function [step, integerstep, ClCd] = GHstep(Cl, Airfoil, airfoil_dir_name)
    %ABCC is Airfoil for which gradient and hessian are being calculated 
    A = Airfoil(1);
    B = Airfoil(2);
    CC = Airfoil(3:4);

    %% Preallocate Arrays
    Airfoils = cell(4);
    clcd = NaN * zeros(4);
    firstpartialderivatives = NaN * zeros(4);
    Gradient = NaN * zeros(3,1);
    Hessian = NaN * zeros(3); %3x3 matrix

    %% Create Arrays
    % Each row contains Airfoils necessary to calculate first partial
    % derivatives at Airfoil in (r,1). Column 1 contains all Airfoils at
    % which first derivatives must be calculated to calculate second 
    % derivatives at Airfoil in (1,1). Row and column 2 represent 
    % derivatives wrt A, row and column 3 wrt B, row and colunm 4 wrt CC.
    
    %  ABCC     (A+1)BCC      A(B+1)CC      AB(CC+1)
    % (A+1)BCC  (A+2)BCC     (A+1)(B+1)CC  (A+1)B(CC+1)
    %  A(B+1)CC (A+1)(B+1)CC  A(B+2)CC      A(B+1)(CC+1)
    %  AB(CC+1) (A+1)B(CC+1)  A(B+1)(CC+1)  AB(CC+2)
    
    for r = 1:4
       if r == 1
           a = A;
           b = B;
           cc = CC;
       elseif r == 2
           a = num2str(str2double(A) + 1);
           b = B;
           cc = CC;
       elseif r == 3
           a = A;
           b = num2str(str2double(B) + 1);
           cc = CC;
       elseif r == 4
           a = A;
           b = B;
           cc = num2str(str2double(CC) + 1);
       end
       
        Airfoils{r,1} = [a b cc];
        Airfoils{r,2} = [num2str(str2double(a)+1) b cc];
        Airfoils{r,3} = [a num2str(str2double(b)+1) cc];
        Airfoils{r,4} = [a b num2str(str2double(cc)+1)];
        
        for c = 1:4
            
            if length(Airfoils{r,c}) == 3 % CC needs to be 2 digits, so put zero in front if missing
                Airfoils{r,c} = [Airfoils{r,c}(1:2), '0', Airfoils{r,c}(3)];
            end
            clcd(r,c) = findclcd2(Cl, Airfoils{r,c},airfoil_dir_name);
            
            %(y2-y1)/(x2-x1), but (x2-x1) always = 1
            firstpartialderivatives(r,c) = clcd(r,c) - clcd(r,1);
            % Comment: there's no point in doing the first column
                % first row is first partial derivatives at ABCC
                % second row is first partial derivatives at (A+1)BCC
                % third row is first partial derivatives at A(B+1)CC
                % fourth row is first partial derivatives at AB(CC+1)
        end
    end
        firstpartialderivatives(:,1) = []; %firt column all zeros
        
        for i = 1:3
            % Gradient is column matrix of all first partial derivatives at
            % Airfoil in (1,1)
            Gradient(i,1) = firstpartialderivatives(1,i);
        end
        for r = 2:4
            for c = 1:3
                % Subtract first row from every other row
                % (y2-y1)/(x2-x1), but (x2-x1) always = 1
                Hessian(r-1,c) = firstpartialderivatives(r,c) - firstpartialderivatives(1,c);
            end
        end
        
        Gradient
        Hessian
        
          ClCd = clcd(1,1); %ClCd at ABCC
          step = Hessian\Gradient %Positive for now because using positive clcd values
          integerstep = fix(step); %Keep A B CC as integers, round step to be smaller
        
%          disp(Airfoils);
%          disp(clcd); 
%          disp(firstpartialderivatives); 
%          disp(Gradient); 
%          disp(Hessian);

end