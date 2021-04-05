function [Cl, best_airfoil, ClCd] = Airfopt2(Cl, GuessAirfoil, airfoil_dir_name)
    %airfoil_dir_name = 'airfoil_database\';
    addpath(airfoil_dir_name);

    if isnumeric(GuessAirfoil)
        trimAirfoil = num2str(GuessAirfoil);
    elseif lower(GuessAirfoil(1:4)) == 'naca'
        trimAirfoil = GuessAirfoil(5:8);
    end
    
    A = str2double(trimAirfoil(1));
    B = str2double(trimAirfoil(2));
    CC = str2double(trimAirfoil(3:end));
    
    if A > 7 || B > 7 || CC > 30
        best_airfoil = 'Improper guess airfoil';
        ClCd = NaN;
    else
        
    G = 1;
    iter = 0;
    maxiter = 20;
    tested = strings(1,maxiter);
    condition = true;
    
        while max(abs(G)) > 0.1 && iter < maxiter && condition
            iter = iter + 1;
            disp(iter); %TEMPORARY FOR TESTING
            Aval(iter) = A;
            Bval(iter) = B;
            CCval(iter) = CC;
            tested(iter) = string([num2str(A), num2str(B), num2str(CC)]);
            disp(tested(iter));%TEMPORARY FOR TESTING
    %         if A == 0 && B ~= 0 %BOUNDARY CASE
    %             B = 0;               
    %         end

%             [G, H, ClCd(iter)] = calc_Gradient_Hessian(Cl, A, B, CC, airfoil_dir_name);
            [G, H, ClCd] = calc_Gradient_Hessian(Cl, A, B, CC, airfoil_dir_name);
            step = -H\G;%-inverseH * G;
            integerstep = fix(step);
%             disp(ClCd(iter));%TEMPORARY FOR TESTING
            disp(step);%TEMPORARY FOR TESTING
            disp(integerstep);%TEMPORARY FOR TESTING
            
            ClCdplot(A,B,CC);

            if A+integerstep(1) >= 0 && A+integerstep(1) <= 7 
                A = A + integerstep(1);
            end
            if B+integerstep(2) >= 0 && B+integerstep(2) <= 7
                B = B + integerstep(2);
            end
            if CC+integerstep(3) >= 0 && CC+integerstep(3) <= 30
                CC = CC + integerstep(3);  
            end
            if iter >= 2 && tested(iter) == tested(iter - 1)
                condition = false;
                tested(tested == "") = [];
            end
            if iter >= 3 && tested(iter) == tested(iter - 2)
                % jumping between 2 airfoils
            end
        best_airfoil = tested(end);
        ClCd = ClCd(end);
        end  
end



