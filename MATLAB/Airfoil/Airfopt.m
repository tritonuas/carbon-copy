function [bestAirfoil, bestClCd] = Airfopt(Cl, GuessAirfoil)

airfoil_dir_name = 'C:\TUAS\carbon-copy\MATLAB\Airfoil\airfoil_database\';

if isnumeric(GuessAirfoil)
    trimAirfoil = num2str(GuessAirfoil);
elseif lower(GuessAirfoil(1:4)) == 'naca'
    trimAirfoil = GuessAirfoil(5:8);
end
A = trimAirfoil(1);
B = trimAirfoil(2);
CC = trimAirfoil(3:end);

% numA = str2double(A);
% numB = str2double(B);
% numCC = str2double(CC);

%N = max([numA, (10 - numA)]);   % Number of iterations cap, only varying A
N = 20;

for p = 1%:N
    if A == '0' && B ~= '0'
    B = '0';               
    fprintf('Using naca%s%s%s because location of max camber must = 0 when max camber = 0\n', A, B, CC);
    end
    [dA(p), ddA(p), ClCdA(p)] = calcderivatives(A, B, CC, Cl, A, airfoil_dir_name);
    [dB(p), ddB(p), ClCdB(p)] = calcderivatives(A, B, CC, Cl, B,  airfoil_dir_name);
    [dCC(p), ddCC(p), ClCdCC(p)] = calcderivatives(A, B, CC, Cl, CC,  airfoil_dir_name);

    gradient(p) = [dA(p); dB(p); dCC(p)];

    hold on
    plot(xA2, yA2, 'r.');
    xlabel('A'); ylabel('Cl/Cd'); 
    title(['Cl/Cd at Cl = ' num2str(Cl) ' vs max camber']);
    axis([0 10 -50 0]);
    
    figure
    plot(xB2, yB2, 'r.');
    xlabel('B'); ylabel('Cl/Cd'); 
    title(['Cl/Cd at Cl = ' num2str(Cl) ' vs location of max camber']);
    axis([0 10 -50 0]);
    
    figure
    plot(xCC2, yCC2, 'r.');
    xlabel('CC'); ylabel('Cl/Cd'); 
    title(['Cl/Cd at Cl = ' num2str(Cl) ' vs max thickness']);
    axis([0 10 -50 0]);
 
end

end

function [derivative, secondderivative, ClCd] = calcderivatives(a, b, cc, CL, v, fdirectory)
    count = 0;
    arange = [a, num2str(str2double(a)+1)];
    brange = [b, num2str(str2double(b)+1)];
    ccrange = [cc, num2str(str2double(cc)+1)];
    if v == 'A' || v == 'a'
        arange = [a, num2str(str2double(a)+1), num2str(str2double(a)+2)];
    elseif v == 'B' || v == 'b'
        brange = [b, num2str(str2double(b)+1), num2str(str2double(b)+2)];
    elseif v == 'CC' || v == 'cc'
        ccrange = [cc, num2str(str2double(cc)+1), num2str(str2double(cc)+2)];
    else
        fprintf('Choose v value of ''A'', ''B'', or ''CC''\n');
    end
    for i = arange 
        for j = brange
            for k = ccrange
            count  = count + 1;
            fID = fopen(strcat(fdirectory, 'naca', i, j, k, '.pol'), 'r');
            if fID == -1 % If file does not exist, generate file
            fID = genXfoil(fdirectory, i, j, k);
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

            data.name{count} = ['naca' i j k];
            if v == 'A' || v == 'a'
                data.variable(count) = str2double(i); %A
            elseif v == 'A' || v == 'a'
                data.variable(count) = str2double(j); %B
            elseif v == 'CC' || v == 'cc'
                data.variable(count) = str2double(k); %CC
            end 
            data.close(count) = Close; % Cl Approximation
            data.clcd(count) = -clcd(cl == Close); % -Cl/Cd at Cl Approximation
            end
        end
    end
   x = data.variable(1);
   xp = data.variable(2);
   xp2 = data.variable(3);
   y = data.clcd(1);
   yp = data.clcd(2);
   yp2 = data.clcd(3);
   
   d2 = (yp2 - yp)/(xp2 - xp);
   
   derivative = (yp - y)/(xp - x); %xp-x = 1
   secondderivative = (d2 - derivative)/(xp - x);
   ClCd = data.clcd(1);
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