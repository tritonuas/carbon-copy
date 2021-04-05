function [Gradient, Hessian, ClCd] = calc_Gradient_Hessian(Cl, A, B, CC, airfoil_dir_name)  
    Airfoils = cell(2,4,3);
    clcd = cell(2,4,3);
    frstprtialderivative = cell(1,4,3);
    Gradient = zeros(3,1);
    Hessian = zeros(3);
    for s = 1:3 % three dimensions for debugging, 3rd index for arrays not necessary
                %1 ==> A; 2 ==> B; 3 ==> CC
        r = 0;
        amove = [A A+1];
        bmove = [B B+1];
        ccmove = [CC CC+1];
        if s == 1
            move = amove;
        elseif s == 2
            move = bmove;
        elseif s == 3
            move = ccmove;
        end
        a = A;
        b = B;
        cc = CC;
        for i = move
            r = r + 1;
            if s == 1
                a = i;
            elseif s == 2
                b = i;
            elseif s == 3
                cc = i;
            end    
            Airfoils{r,1,s} = [num2str(a), num2str(b),num2str(cc)];
            Airfoils{r,2,s} = [num2str(a+1), num2str(b), num2str(cc)];
            Airfoils{r,3,s} = [num2str(a), num2str(b+1), num2str(cc)];
            Airfoils{r,4,s} = [num2str(a), num2str(b), num2str(cc+1)];
%                  disp(Airfoils); %TEMPORARY FOR TESTING
            for c = 1:4
                clcd{r,c,s} = findclcd(Cl, Airfoils{r,c,s}, airfoil_dir_name);
%                 disp(clcd); %TEMPORARY FOR TESTING
            end
        end
        for c = 1:4
            frstprtialderivative{1,c,s} = (-clcd{2,c,s}) - (-clcd{1,c,s});% negatives to make optimization find minimum
        end
        
        Gradient(s,1) = frstprtialderivative{1,1,s};
        for clm = 2:4
            Hessian(s,clm-1) = frstprtialderivative{1,clm,s} - frstprtialderivative{1,1,s}; %Hessian needs some work
        end
    end
    ClCd = clcd{1,1,1};
%      disp(Gradient);disp(Hessian);%TEMPORARY FOR TESTING
end