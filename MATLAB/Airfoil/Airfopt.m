function bestAirfoil = Airfopt(Cl, GuessAirfoil)

airfoil_dir_name = 'C:\TUAS\carbon-copy\MATLAB\Airfoil\airfoil_database\';

if isnumeric(GuessAirfoil)
    trimAirfoil = num2str(GuessAirfoil);
elseif lower(GuessAirfoil(1:4)) == 'naca'
    trimAirfoil = GuessAirfoil(5:8);
end
A = trimAirfoil(1);
B = trimAirfoil(2);
CC = trimAirfoil(3:end);
if A == '0' && B ~= '0'
    B = '0';               
    fprintf('Using naca%s%s%s because location of max camber must = 0 when max camber = 0\n', A, B, CC);
end
numA = str2double(A);
% numB = str2double(B);
% numCC = str2double(CC);

N = max([numA, (10 - numA)]);   % Number of iterations cap, only varying A
data = struct('Airfoil', cell(N,2), 'Cl', cell(N,2), 'Close', cell(N,2), 'A', cell(N,2), 'ClCd', 0);
for i = 1:N
    count = 0;
    for j = [A, num2str(numA + 1)]
        count  = count + 1;
        fID = fopen(strcat(airfoil_dir_name, 'naca', j, B, CC, '.pol'), 'r');
        if fID == -1 % If file does not exist, generate file
            fID = genXfoil(airfoil_dir_name, j, B, CC);
        end
        D = textscan(fID,'%f %f %f %f %f %f %f', 'HeaderLines', 12);  %skips headers
        fclose(fID);
        
%         alpha = D{:,1};                     %Alpha
        cl = D{:,2};                        %Coefficient of Lift
        cd = D{:,3};                        %Coefficient of Drag
        clcd = cl./cd;                      %drag efficiency

        error = abs(cl-Cl);                 % What if there are two places where cl == Cl
        Close = cl(error == min(error));
        Close = Close(1);

        fclose('all');
        
        
        data(i).Airfoil = ['naca' j B CC];
        data(i).A(count) = str2double(j);
        data(i).Cl = Cl;
        data(i).Close(count) = Close;
        data(i).ClCd(count) = clcd(cl == Close);
    end

    x1 = data(i).A(1);
    x2 = data(i).A(2);
    y1 = data(i).ClCd(1);
    y2 = data(i).ClCd(2);
    dA = (y2 - y1)/ (x2 - x1); % x2-x1 always = 1
    
    hold on
    plot(x1, y1, 'r.');
    xlabel('A'); ylabel('Cl/Cd'); 
    title(['Cl/Cd at Cl = ' num2str(Cl) ' vs max camber']);
    axis([0 10 0 50]);

    if dA > 0 && numA < 9 
        numA = numA + 1; 
    elseif dA < 0 && numA > 0 
        numA = numA - 1;
    end
 
 A = num2str(numA);
end
bestAirfoil = data.Airfoil(data.ClCd == max(data.ClCd(1)));

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