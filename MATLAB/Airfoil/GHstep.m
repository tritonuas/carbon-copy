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
%     for s = 1:3 % three dimensions for debugging, 3rd index for arrays not necessary
%                 %1 ==> A; 2 ==> B; 3 ==> CC
%         r = 0;
%         amove = [A A+1];
%         bmove = [B B+1];
%         ccmove = [CC CC+1];
%         if s == 1
%             move = amove;
%         elseif s == 2
%             move = bmove;
%         elseif s == 3
%             move = ccmove;
%         end
%         a = A;
%         b = B;
%         cc = CC;
%         for i = move
%             r = r + 1;
%             if s == 1
%                 a = i;
%             elseif s == 2
%                 b = i;
%             elseif s == 3
%                 cc = i;
%             end    
%             Airfoils{r,1} = [num2str(a), num2str(b),num2str(cc)];
%             Airfoils{r,2} = [num2str(a+1), num2str(b), num2str(cc)];
%             Airfoils{r,3} = [num2str(a), num2str(b+1), num2str(cc)];
%             Airfoils{r,4} = [num2str(a), num2str(b), num2str(cc+1)];
% %                  disp(Airfoils); %TEMPORARY FOR TESTING
%             for c = 1:4
%                 clcd{r,c} = findclcd2(Cl, Airfoils{r,c}, airfoil_dir_name);
% %                 disp(clcd); %TEMPORARY FOR TESTING
%             end
%         end
%         for c = 1:4
%             frstprtialderivative{1,c} = (-clcd{2,c}) - (-clcd{1,c});% negatives to make optimization find minimum
%         end
%         
%         Gradient(s,1) = frstprtialderivative{1,1};
%         for clm = 2:4
%             Hessian(s,clm-1) = frstprtialderivative{1,clm} - frstprtialderivative{1,1}; %Hessian needs some work
%         end
%     end
    
    %% Brute Force Airfoil Array, doesn't work now that A, B, CC are strings
%     Airfoils{1,1} = [num2str(A) num2str(B) num2str(CC)];
%     Airfoils{1,2} = [num2str(A+1) num2str(B) num2str(CC)];
%     Airfoils{1,3} = [num2str(A) num2str(B+1) num2str(CC)];
%     Airfoils{1,4} = [num2str(A) num2str(B) num2str(CC+1)];
%     
%     Airfoils{2,1} = [num2str(A+1) num2str(B) num2str(CC)];
%     Airfoils{2,2} = [num2str(A+2) num2str(B) num2str(CC)];
%     Airfoils{2,3} = [num2str(A+1) num2str(B+1) num2str(CC)];
%     Airfoils{2,4} = [num2str(A+1) num2str(B) num2str(CC+1)];
%     
%     Airfoils{3,1} = [num2str(A) num2str(B+1) num2str(CC)];
%     Airfoils{3,2} = [num2str(A+1) num2str(B+1) num2str(CC)];
%     Airfoils{3,3} = [num2str(A) num2str(B+2) num2str(CC)];
%     Airfoils{3,4} = [num2str(A) num2str(B+1) num2str(CC+1)];
%     
%     Airfoils{4,1} = [num2str(A) num2str(B) num2str(CC+1)];
%     Airfoils{4,2} = [num2str(A+1) num2str(B) num2str(CC+1)];
%     Airfoils{4,3} = [num2str(A) num2str(B+1) num2str(CC+1)];
%     Airfoils{4,4} = [num2str(A) num2str(B) num2str(CC+2)];
    
    %% More elegant Airfoil Array
    for r = 1:4
       if r == 1                    %Some conversion between str and double unnecesary, done for consistency
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
       
       % Each row has airfoils necessary to calculate first partial derivative at Airfoil (r,1) 
        Airfoils{r,1} = [a b cc];
        Airfoils{r,2} = [num2str(str2double(a)+1) b cc];
        Airfoils{r,3} = [a num2str(str2double(b)+1) cc];
        Airfoils{r,4} = [a b num2str(str2double(cc)+1)];

        for c = 1:4
            
            if length(Airfoils{r,c}) == 3 % CC needs to be 2 digits, str2double erases leading zero
                Airfoils{r,c} = [Airfoils{r,c}(1:2), '0', Airfoils{r,c}(3)];
            end
            
            clcd(r,c) = findclcd2(Cl, Airfoils{r,c},'airfoil_database\\specific_cl\\');
            firstpartialderivatives(r,c) = clcd(r,c) - clcd(r,1); 
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
        
        disp(Airfoils);
        disp(clcd); 
%         disp(firstpartialderivatives); 
%         disp(Gradient); 
        disp(Hessian);

end