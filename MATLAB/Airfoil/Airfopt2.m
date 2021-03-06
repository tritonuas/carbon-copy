function [bestAirfoil, bestClCd] = Airfopt2(Cl, GuessAirfoil)
    airfoil_dir_name = 'C:\TUAS\carbon-copy\MATLAB\Airfoil\airfoil_database\';

    if isnumeric(GuessAirfoil)
        trimAirfoil = num2str(GuessAirfoil);
    elseif lower(GuessAirfoil(1:4)) == 'naca'
        trimAirfoil = GuessAirfoil(5:8);
    end
    A = trimAirfoil(1);
    B = trimAirfoil(2);
    CC = trimAirfoil(3:end);
    
    for p = 1% Will eventually be a while loop referencing gradient
        if A == '0' && B ~= '0'
            B = '0';               
            fprintf('Using naca%s%s%s because location of max camber must = 0 when max camber = 0\n', A, B, CC);
        end
        Airfoils = cell(2,4,3);
        negclcd = cell(2,4,3);
        pder = cell(1,4,3);
        for v = 1:3
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
                Airfoils{r,1,v} = [num2str(a), num2str(b),num2str(cc)];
                Airfoils{r,2,v} = [num2str(a+1), num2str(b), num2str(cc)];
                Airfoils{r,3,v} = [num2str(a), num2str(b+1), num2str(cc)];
                Airfoils{r,4,v} = [num2str(a), num2str(b), num2str(cc+1)];

                for c = 1:4
                    negclcd{r,c,v} = findclcd(Cl, Airfoils{r,c,v}, airfoil_dir_name);
                end
            end
            for c = 1:4
                pder{1,c,v} = negclcd{2,c,v} - negclcd{1,c,v};
                
            end
            grad(v,1) = pder{1,1,v};
        end
    end
end

function nClCd = findclcd(CL, Airfoil, fdirectory)
    fID = fopen(strcat(fdirectory, 'naca', Airfoil, '.pol'), 'r');
    if fID == -1 % If file does not exist, generate file
        fID = genXfoil(fdirectory, Airfoil);
    end
    D = textscan(fID,'%f %f %f %f %f %f %f', 'HeaderLines', 12);  %skips headers
    fclose(fID);

    %         alpha = D{:,1};                     %Alpha
    cl = D{:,2};                        %Coefficient of Lift
    cd = D{:,3};                        %Coefficient of Drag
    clcd = cl./cd;                      %drag efficiency

    error = abs(cl-CL);                 % What if there are two places where cl == Cl
    Close = cl(error == min(error));
    Close = Close(1);

    fclose('all');

    nClCd = -clcd(cl == Close); % -Cl/Cd at Cl Approximation
end

function FID = genXfoil(directory, A, B, CC)
    airfoil_d_name = "airfoil_database";
    addpath(airfoil_d_name);
    strAirfoil = [A B CC];

    fid = fopen('xfoil_input.txt','w');   % create the inputs  
    fprintf(fid,["naca "+strAirfoil+"\n"]);
    fprintf(fid,'pane\n');   % makes the airfoils nice
    fprintf(fid,['oper\n']);fprintf(fid,['visc\n']);   % makes visc analysis w/ Re=4E5
    fprintf(fid,['4e5\n']);fprintf(fid,['pacc\n']);
    fprintf(fid,"airfoil_database/naca"+strAirfoil+".pol\n\n");
    fprintf(fid,'iter\n');fprintf(fid,'25\n'); % aint nobody got time fo 50
    fprintf(fid,'aseq 0 20 0.5\n');% alpha 0 to 20 deg in 0.5 deg increment
    fprintf(fid,'pacc\n');

    cmd = 'xfoil.exe < xfoil_input.txt';   % running on xfoil
    [status,result] = system(cmd);

    fclose('all');
    delete('xfoil_input.txt');
    FID = fopen(strcat(directory, 'naca', strAirfoil, '.pol'), 'r');
end