function [step, ClCd] = simplestep(Cl, Airfoil, airfoil_dir_name)

    A = Airfoil(1);
    B = Airfoil(2);
    CC = Airfoil(3:4);

    %% Preallocate Arrays
    Airfoils = cell(1,4);
    clcd = NaN * zeros(1,4);
    firstpartialderivatives = NaN * zeros(1,4);
    Gradient = NaN * zeros(3,1);
    step = NaN * zeros(3,1);
    
     %% Create Arrays
        
        Airfoils{1,1} = [A B CC];
        Airfoils{1,2} = [num2str(str2double(A)+1) B CC];
        Airfoils{1,3} = [A num2str(str2double(B)+1) CC];
        Airfoils{1,4} = [A B num2str(str2double(CC)+1)];

        for c = 1:4
            
            if length(Airfoils{1,c}) == 3 % CC needs to be 2 digits, str2double erases leading zero
                Airfoils{1,c} = [Airfoils{1,c}(1:2), '0', Airfoils{1,c}(3)];
            end
            
            clcd(1,c) = findclcd2(Cl, Airfoils{1,c},'airfoil_database\\specific_cl\\');
            firstpartialderivatives(1,c) = clcd(1,c) - clcd(1,1);
        end
        firstpartialderivatives(1) = [];
        
        for i = 1:3
            Gradient(i,1) = firstpartialderivatives(1,i);
            if Gradient(i) > 0
                step(i) = 1;
            elseif Gradient(i) < 0
                step(i) = -1;
            else
                step(i) = 0;
            end
        end
        
        ClCd = clcd(1);
        
%          disp(Airfoils);
%          disp(clcd); 
%          disp(firstpartialderivatives); 
%          disp(Gradient); 
end
        