function [step] = Airfopt2(Cl, GuessAirfoil)
    airfoil_dir_name = 'airfoil_database\';

    if isnumeric(GuessAirfoil)
        trimAirfoil = num2str(GuessAirfoil);
    elseif lower(GuessAirfoil(1:4)) == 'naca'
        trimAirfoil = GuessAirfoil(5:8);
    end
    
    A = str2double(trimAirfoil(1));
    B = str2double(trimAirfoil(2));
    CC = str2double(trimAirfoil(3:end));
    
    G = 1;
    iter = 0;
    maxiter = 20;
    tested = strings(1,maxiter);
    
    while max(abs(G)) > 0.1 && iter < maxiter
        iter = iter + 1;
        tested(iter) = string([num2str(A), num2str(B), num2str(CC)]);
        
        if A == 0 && B ~= 0
            B = 0;               
            fprintf('Using naca%s%s%s because location of max camber must = 0 when max camber = 0\n', A, B, CC);
        end
        [G, H] = calc_Gradient_Hessian(Cl, A, B, CC, airfoil_dir_name);
        
        step = H\-G;
        integerstep = fix(step);
        
        if A+integerstep(1) >= 0 && A+integerstep(1) <= 8 
        A = A + integerstep(1);
        end
        if B+integerstep(2) >= 0 && B+integerstep(2) <= 8
        B = B + integerstep(2);
        end
        if CC+integerstep(3) >= 0 && CC+integerstep(3) <= 30
        CC = CC + integerstep(3);  
        end
    end
    
    final_airfoil = tested;
    disp(iter);
    disp("final_airfoil: " + final_airfoil);
    
end



