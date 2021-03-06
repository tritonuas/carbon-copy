function [Gradient, Hessian] = calc_Gradient_Hessian2(Cl, A, B, CC, directory)  
    Airfoils = cell(4,4);
    negclcd = cell(4,4);
    frstprtialderderivative = cell(3,4);
    Gradient = zeros(3,1);
    Hessian = zeros(3);
    Airfoils{1,1} = [num2str(A), num2str(B),num2str(CC)];
    for v = 1:3 % three dimensions for debugging, 3rd index for arrays not necessary
                %1 ==> A; 2 ==> B; 3 ==> CC
        r = 0;
        amove = [A A+1];
        bmove = [B B+1];
        ccmove = [CC CC+1];
        if v == 1
            move = amove;
        elseif v == 2
            move = bmove;
        elseif v == 3
            move = ccmove;
        end
        a = A;
        b = B;
        cc = CC;
        for i = move
            r = r + 1;
            if v == 1
                a = i;
            elseif v == 2
                b = i;
            elseif v == 3
                cc = i;
            end    
            
            Airfoils{r,2} = [num2str(a+1), num2str(b), num2str(cc)];
            Airfoils{r,3} = [num2str(a), num2str(b+1), num2str(cc)];
            Airfoils{r,4} = [num2str(a), num2str(b), num2str(cc+1)];
                 disp(Airfoils); %TEMPORARY FOR TESTING
            for c = 1:4
                negclcd{r,c} = findclcd(Cl, Airfoils{r,c}, directory);
            end
        end
        for c = 1:4
            frstprtialderderivative{1,c} = negclcd{2,c} - negclcd{1,c};
        end
        
        Gradient(v,1) = frstprtialderderivative{1,1};
        for clm = 2:4
            Hessian(v,clm-1) = frstprtialderderivative{1,clm} - frstprtialderderivative{1,1}; %Hessian needs some work
        end
    end
     disp(Gradient);disp(Hessian);%TEMPORARY FOR TESTING
end