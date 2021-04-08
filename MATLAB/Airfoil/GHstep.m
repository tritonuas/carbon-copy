function [step, integerstep, ClCd] = GHstep(Cl, Airfoil, airfoil_dir_name)
    
    A = Airfoil(1);
    B = Airfoil(2);
    CC = Airfoil(3:4);

    %% Preallocate Arrays
    Airfoils = cell(2,4);
    clcd = NaN * zeros(2,4);
    firstpartialderivatives = NaN * zeros(3,4);
    Gradient = NaN * zeros(3,1);
    Hessian = NaN * zeros(3);

    %% Create Arrays
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
       
       % Each row has airfoils necessary to calculate first partial derivatives at Airfoil (r,1) 
        Airfoils{r,1} = [a b cc];
        Airfoils{r,2} = [num2str(str2double(a)+1) b cc];
        Airfoils{r,3} = [a num2str(str2double(b)+1) cc];
        Airfoils{r,4} = [a b num2str(str2double(cc)+1)];

        for c = 1:4
            
            if length(Airfoils{r,c}) == 3 % CC needs to be 2 digits, str2double erases leading zero
                Airfoils{r,c} = [Airfoils{r,c}(1:2), '0', Airfoils{r,c}(3)];
            end
            
            clcd(r,c) = findclcd2(Cl, Airfoils{r,c},'airfoil_database\\specific_cl\\');
            firstpartialderivatives(r,c) = (-clcd(r,c)) - (-clcd(r,1)); % negative so optimization minimizes
        end
    end
        firstpartialderivatives(:,1) = [];
        
        for i = 1:3
            Gradient(i,1) = firstpartialderivatives(1,i);
        end
        for r = 2:4
            for c = 1:3
                Hessian(r-1,c) = firstpartialderivatives(r,c) - firstpartialderivatives(1,c);
            end
        end
        
          ClCd = clcd(1,1);
          step = -Hessian\Gradient;
          integerstep = fix(step);
        
%         disp(Airfoils);
%         disp(clcd); 
%         disp(firstpartialderivatives); 
%         disp(Gradient); 
%         disp(Hessian);

end