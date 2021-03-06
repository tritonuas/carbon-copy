function [bestAirfoil, bestClCd] = Airfopt2(Cl, GuessAirfoil)
    airfoil_dir_name = 'C:\TUAS\carbon-copy\MATLAB\Airfoil\airfoil_database\';

    if isnumeric(GuessAirfoil)
        trimAirfoil = num2str(GuessAirfoil);
    elseif lower(GuessAirfoil(1:4)) == 'naca'
        trimAirfoil = GuessAirfoil(5:8);
    end
    A = str2double(trimAirfoil(1));
    B = str2double(trimAirfoil(2));
    CC = str2double(trimAirfoil(3:end));
    
    for p = 1% Will eventually be a while loop referencing gradient
        if A == 0 && B ~= 0
            B = 0;               
            fprintf('Using naca%s%s%s because location of max camber must = 0 when max camber = 0\n', A, B, CC);
        end
        [G, H] = calc_Gradient_Hessian(Cl, A, B, CC, airfoil_dir_name);
        
        step = -G./H;
        
    end
end



