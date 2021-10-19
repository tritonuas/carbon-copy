%% Airfoil Optimization
% curly brackets so dims dont gotta match up
% Erik Lovekin 1/2021
% Graphing/File Access/File Generation based on Kaelen Tan's Airfoil Comparison 2
% Can't handle files that are not in a designation+number naming format
% Only changes max camber
% Uses slope between two points instead of gradient or derivative at one point
% Functions to look in to: strtok, isstrprop

% Make use of cell arrays instead of converting between double and string
%% Function Name
function best_airfoil = Optfoil(Cl, GuessAirfoil) %, clStep)
%% General
close all;
format loose;
format shortg;

%% File directory stuff
airfoil_dir_name = 'airfoil_database/'; 
addpath(airfoil_dir_name); % Don't think I'm using this right

%% Define Starting Conditions
% Cl = 0.9673;      %Improvement: make it able to accept either 'naca2212' or '2212'
% trimAirfoil = GuessAirfoil;
A = GuessAirfoil(5);
numA = str2double(A);
CC = GuessAirfoil(7:end);
B = GuessAirfoil(6);
if A == '0' && B ~= '0'
    B = '0';                    % Only allow location of max camber to be zero when max camber = 0
    fprintf('Using naca%s%s%s because location of max camber must = 0 when max camber = 0\n', A, B, CC);
end

N = max([numA, (10 - numA)]);   % Number of iterations cap, only varying A
E = cell(N, 3);
i = 1;

% P = zeros(N,4) * NaN;
% slope = 2;
%% Break needs to leave
    for j = 1:N
        strAirfoil = strcat(A, B, CC);
        numA = str2double(A);
        numB = str2double(B);
        numCC = str2double(CC);
              
%% Reading data from file

        fID = fopen(strcat(airfoil_dir_name, 'naca', strAirfoil, '.pol'), 'r');
        if fID == -1 % If file does not exist, generate file

            % Set up the file names
            airfoil_d_name = "airfoil_database";
            addpath(airfoil_d_name);

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
            fID = fopen(strcat(airfoil_dir_name, 'naca', strAirfoil, '.pol'), 'r');                
        end

        D = textscan(fID,'%f %f %f %f %f %f %f', 'HeaderLines', 12);  %skips headers
        fclose(fID);

        alpha = D{:,1};                     %Alpha
        cl = D{:,2};                        %Coefficient of Lift
        cd = D{:,3};                        %Coefficient of Drag
        clcd = cl./cd;                      %drag efficiency

        %                     test = [alpha, cl, clcd];
        %                     disp(test);

        error = abs(cl-Cl);                 % What if there are two places where cl == Cl
        Close = cl(error == min(error));
        Close = Close(1);

        fclose('all');

        aStep = 0.5; %a(max(alpha) - min(alpha))/ (length(alpha)-1);                                            

        alpha_fill = alpha;                                        %Use data directly from file
        cl_fill = cl;
        clcd_fill = clcd; 

        hold on;
        plot(alpha, clcd, 'r.-');
        xlabel('alpha'); ylabel('Cl/Cd'); title('Cl/Cd vs. Alpha');
        plot(alpha(cl == Close), clcd(cl == Close), 'b.', 'Markersize', 20)

        %% Calculate Slope
        alpha_p = alpha_fill(cl_fill == Close);
        clcd_p = clcd_fill(cl_fill == Close);
        x1 = alpha_p;
        x2 = alpha_p + aStep;
        y1 = clcd_p;
        y2 = clcd_fill(alpha_fill == x2(1)); % naca2118
        slope = (y2 - y1(1))/(x2(1) - x1(1));

        E(j,:) = {strcat('naca',strAirfoil), slope, clcd_p};

        if slope < 0 && numA < 9
            i = 1;
        elseif slope > 0 && numA > 0
            i = -1;
        else
            i = 0;
        end

        numA = numA + i;

        A = num2str(numA);
        B = num2str(numB);
        CC = num2str(numCC);   
    end
    
    emptyCells = cellfun('isempty', E); 
    E(all(emptyCells,2),:) = [];
    disp(E);

    [mx, in] = max([E{:,3}]);
    best_airfoil = E{in, 1};
    disp(best_airfoil);
    disp(mx);
        
%% Plotting

fID = fopen(strcat(airfoil_dir_name, best_airfoil, '.pol'), 'r');
D = textscan(fID,'%f %f %f %f %f %f %f', 'HeaderLines', 12);  %skips headers
fclose(fID);

alpha = D{:,1};                     %Alpha
cl = D{:,2};                        %Coefficient of Lift
cd = D{:,3};                        %Coefficient of Drag
clcd = cl./cd;                      %drag efficiency

%                     test = [alpha, cl, clcd];
%                     disp(test);

error = abs(cl-Cl);                 % What if there are two places where cl == Cl
Close = cl(error == min(error));
Close = Close(1);

fclose('all');
plot(alpha, clcd, 'k.-');
plot(alpha(cl == Close), clcd(cl == Close), 'bo', 'Markersize', 10, 'Linewidth', 3)
end