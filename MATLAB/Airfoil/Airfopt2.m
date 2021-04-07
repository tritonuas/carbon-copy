
function [best_airfoil, ClCd] = Airfopt2(Cl, GuessAirfoil, airfoil_dir_name)
%   airfoil_dir_name = 'airfoil_database\\specific_cl\\';
%   GuessAirfoil of format 'nacaABCC';
    addpath(airfoil_dir_name);

    %% Standardize Input
%     if isnumeric(GuessAirfoil)
%         trimAirfoil = num2str(GuessAirfoil);
%     elseif lower(GuessAirfoil(1:4)) == 'naca'
%         trimAirfoil = GuessAirfoil(5:8);
%     end

    A = GuessAirfoil(5);
    B = GuessAirfoil(6);
    CC = GuessAirfoil(7:8);
    maxA = 7;
    maxB = 7;
    maxCC = 30;
    minCC = 10;

    if str2double(A) > maxA || str2double(B) > maxB || str2double(CC) > maxCC || str2double(CC) < minCC
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
%             if length(CC) == 1 %Not necessary with minCC >= 10
%                 CC = ['0' CC];
%             end
            
            iter = iter + 1;
%             disp(iter); %TEMPORARY FOR TESTING

            %% For Graphing
            aSet(iter) = str2double(A);
            bSet(iter) = str2double(B);
            ccSet(iter) = str2double(CC);  %OK if doesn't have leading zero
            %%
            
            tested(iter) = [A B CC];
%             disp(tested(iter));%TEMPORARY FOR TESTING

            [~, integerstep, ClCd] = GHstep(Cl, tested(iter), airfoil_dir_name);
            
            ClCdplot(iter,aSet,bSet,ccSet,ClCd,maxA,maxB,maxCC);

            %% Update Values for Next Iteration
            if str2double(A)+integerstep(1) >= 0 && str2double(A)+integerstep(1) <= maxA 
                A = num2str(str2double(A) + integerstep(1));
            end
            if str2double(B)+integerstep(2) >= 0 && str2double(B)+integerstep(2) <= maxB
                B = num2str(str2double(B) + integerstep(2));
            end
            if str2double(CC)+integerstep(3) >= minCC && str2double(CC)+integerstep(3) <= maxCC
                CC = num2str(str2double(CC) + integerstep(3));  
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



