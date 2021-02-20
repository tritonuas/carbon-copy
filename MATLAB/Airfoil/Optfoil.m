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
%  clear all; clc;
format loose;
format shortg;

%% File directory stuff
airfoil_dir_name = 'C:\TUAS\carbon-copy\MATLAB\Airfoil\airfoil_database\'; 
addpath(airfoil_dir_name);

%% Define Starting Conditions
% Cl = 0.9673;      %Improvement: make it able to accept either 'naca2212' or '2212'
% trimAirfoil = GuessAirfoil;
A = GuessAirfoil(5);
initialA = A;
numA = str2double(A);
CC = GuessAirfoil(7:end);
if numA == 0
    B = '0';                    % Only allow location of max camber to be zero when max camber = 0
    fprintf('Using naca%s%s%s because location of max camber must = 0 when max camber = 0\n', A, B, CC);
else
    B = GuessAirfoil(6);
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
              if isfile(strcat(airfoil_dir_name, 'naca', strAirfoil, '.pol'))
                    fID = fopen(strcat(airfoil_dir_name, 'naca', strAirfoil, '.pol'), 'r');%concatenates guess_airfoil and .pol, opens that file for reading
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
                    
%                      cl_range = max(cl) - min(cl);
%                      alpha_range = max(alpha) - min(cl);
%                      clcd_range = max(clcd) - min(clcd);
                        
%% Allow input of any Cl Value (might not be necessary anymore)
%                         if nargin == 3                                                %If clStep is entered
%                             aStep = (cl_range / alpha_range) * clStep;              %Resolution of alpha
%                             clcdStep = (cl_range/ clcd_range) * clStep;             %Resolution of Cl/Cd
%                             
%                             alpha_fill = min(alpha) : aStep : max(alpha);           %Alpha filled in data points
%                             
%                             trans1 = clStep * 10 .^ [1:10];                         %Calculates how many decimal places to round to
%                             rnd = find(trans1 == round(trans1), 1);
%                             cl_fill = round(interp1(alpha, cl, alpha_fill), rnd);   %Cl filled in data points
%                             
%                 ERROR HERE: clcdSTEP is calculated, so not
%                 finite number of decimal places, as was the
%                 case for cl.
%                             trans2 = clcdStep * 10 .^ [1:10]; 
%                             rnd2 = find(trans2 == round(trans2), 1);
%                             clcd_fill = round(interp1(alpha, clcdEfficiency, alpha_fill), rnd);% Cl/Cd filled in data points
%                         else
%% If clStep is not entered, use data directly from file
                            aStep = 0.5; %a(max(alpha) - min(alpha))/ (length(alpha)-1);                                            
                            
                            alpha_fill = alpha;                                        %Use data directly from file
                            cl_fill = cl;
                            clcd_fill = clcd;
%                         end   

%% Calculate Slope
                     alpha_p = alpha_fill(cl_fill == Close);
                     clcd_p = clcd_fill(cl_fill == Close);
                     x1 = alpha_p;
                     x2 = alpha_p + aStep;
                     y1 = clcd_p;
                     y2 = clcd_fill(alpha_fill == x2);
                     slope = (y2 - y1)/(x2 - x1);
                     
                     E(j,:) = {strcat('naca',strAirfoil), slope, clcd_p};
%                      P(j,:) = [strcat('naca',strAirfoil), alpha, cl,clcd]; % For plot at the end 
                     
                 if slope < -1 && numA < 9               % Make slope cuttoff less arbitrary (second derivative?)
                    i = 1;
                 elseif slope > 1 && numA > 0
                    i = -1;
                 else
                     i = 0;
                 end
                 
                 numA = numA + i;
                 
              else          % If file does not exist, generate file
                %%
                % Set up the file names
                airfoil_dir_name = "airfoil_database";
                addpath(airfoil_dir_name);
               
                num = strAirfoil;

                fid = fopen('xfoil_input.txt','w');   % create the inputs  
                fprintf(fid,["naca "+num+"\n"]);
                fprintf(fid,'pane\n');   % makes the airfoils nice
                fprintf(fid,['oper\n']);fprintf(fid,['visc\n']);   % makes visc analysis w/ Re=4E5
                fprintf(fid,['4e5\n']);fprintf(fid,['pacc\n']);
                fprintf(fid,"airfoil_database/naca"+num+".pol\n\n");
                fprintf(fid,'iter\n');fprintf(fid,'25\n'); % aint nobody got time fo 50
                fprintf(fid,'aseq 0 20 0.5\n');% alpha 0 to 20 deg in 0.5 deg increment
                fprintf(fid,'pacc\n');

                cmd = 'xfoil.exe < xfoil_input.txt';   % running on xfoil
                [status,result] = system(cmd);

                fclose('all');
                delete('xfoil_input.txt');
             end
                 
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
% figure;
% hold on;
% plot(P(:,2),P(:,4),'k-');
% xlabel('Alpha'); ylabel('Cl/Cd'); title('Cl/Cd vs. Alpha');
end

%     % CL vs Alpha for data
% %     figure;
% %     hold on;
% %     plot(alpha,cl, '.-', 'MarkerSize', 10);
% %     xlabel('alpha'); ylabel('CL'); title('CL vs. Alpha');
% %     plot(alpha_pass, Cl, '.', 'MarkerSize', 15);
% % 
% %     % CL/CD vs Alpha for data
% %     figure;
% %     hold on
% %     plot(alpha, clcdEfficiency, 'b.-', 'MarkerSize', 10);
% %     xlabel('alpha');ylabel('CL/CD');title('CL/CD vs. alpha');
% %     legend(GuessAirfoil);
% %     plot(alpha_pass, clcdEfficiency(alpha == alpha_pass), 'r.', 'MarkerSize', 15)
% %     plot(alpha_pass+.5,clcdEfficiency(alpha == alpha_pass+.5), 'g.', 'MarkerSize', 15);
% %     plot(alpha_pass-.5,clcdEfficiency(alpha == alpha_pass-.5), 'g.', 'MarkerSize', 15);

% %             cl_fill = min(cl):.0001:max(cl);
% %             alpha_fill = interp1(cl, alpha, cl_fill);
% %             B = [alpha_fill; cl_fill]';
% %             alpha_f = B(:,1);
% %             cl_f = B(:,2);
% %             alpha_p = alpha_f(cl_f == Cl);
% 
%     % CL vs Alpha Interpolated
% %     figure;
% %     hold on;
% %     plot(alpha_f,cl_f, '.-', 'MarkerSize', 10);
% %     xlabel('alpha'); ylabel('CL'); title('CL vs. Alpha');
% %     plot(alpha_p, Cl, '.', 'MarkerSize', 15);

