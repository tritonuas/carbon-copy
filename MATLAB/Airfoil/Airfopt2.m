
function [best_airfoil, ClCd] = Airfopt2(Cl, GuessAirfoil, airfoil_dir_name)
%   airfoil_dir_name = 'airfoil_database\\specific_cl\\';
    addpath(airfoil_dir_name);

    %% Standardize Input
    if isnumeric(GuessAirfoil)
        trimAirfoil = num2str(GuessAirfoil);
    elseif lower(GuessAirfoil(1:4)) == 'naca'
        trimAirfoil = GuessAirfoil(5:8);
    end
    
    A = str2double(trimAirfoil(1));
    B = str2double(trimAirfoil(2));
    CC = str2double(trimAirfoil(3:end));
    maxA = 7;
    maxB = 7;
    maxCC = 30;
    
    if A > maxA || B > maxB || CC > maxCC
        best_airfoil = 'Improper guess airfoil';
        ClCd = NaN;
    else
        %% Prealocating Arrays
        G = 1;
        iter = 0;
        maxiter = 20;
        aSet = NaN * zeros(maxiter,1);
        bSet = NaN * zeros(maxiter,1);
        ccSet = NaN * zeros(maxiter,1);
        tested = strings(maxiter,1);
        ClCd = NaN * zeros(maxiter,1);
        bounce = false;
        
        while max(abs(G)) > 0.1 && iter < maxiter && ~bounce
            iter = iter + 1;
%             disp(iter); %TEMPORARY FOR TESTING
            aSet(iter) = A;
            bSet(iter) = B;
            ccSet(iter) = CC;
            tested(iter) = string([num2str(A), num2str(B), num2str(CC)]);
%             disp(tested(iter));%TEMPORARY FOR TESTING

            [~, integerstep] = GHstep(Cl, A, B, CC, airfoil_dir_name);
            
%             ClCdplot(iter,aSet,bSet,ccSet,ClCd,maxA,maxB,maxCC);

            if A+integerstep(1) >= 0 && A+integerstep(1) <= 7 
                A = A + integerstep(1);
            end
            if B+integerstep(2) >= 0 && B+integerstep(2) <= 7
                B = B + integerstep(2);
            end
            if CC+integerstep(3) >= 0 && CC+integerstep(3) <= 30
                CC = CC + integerstep(3);  
            end
            %% Test for Bouncing Between Airfoils
%             if iter >= 2 && tested(iter) == tested(iter - 1) % stuck on one airfoil
%                 bounce = true;
%                 tested(tested == "") = [];
%                 ClCd(isnan(ClCd)) = [];
%             end
%             if iter >= 3 && tested(iter) == tested(iter - 2) % jumping between 2 airfoils
%                 bounce = true;
%                 tested(tested == "") = [];
%                 ClCd(isnan(ClCd)) = [];
%             end
        end
        hold off
        disp(tested);
        best_airfoil = tested(end);
    end  

end



