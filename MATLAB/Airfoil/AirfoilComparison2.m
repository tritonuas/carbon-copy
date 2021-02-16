%% Airfoil Comparison 2
% curly brackets so dims dont gotta match up
% Kaelan Tan 11/2020
% Based on Andrew Fletcher's Airfoil Comparison
% Edits by Andrew Fletcher
% Can't really handle files that are not in a designation+number naming
% format

%% General
clear,clc;
close all;
clc;
format compact;
format shortG;

%% File directory stuff
A=dir('airfoil_database');      % your folder of airfoils
addpath("airfoil_database");
names={A.name}; names(1)=[]; names(1)=[];

%% Inputs
cont=0;
i=1;
airfoil_name = "initial_string";
while (cont~='N' && cont~='n')
    airfoil_name(i) = input("Airfoil Name: ",'s');
%     num(i)=input("Airfoil Number: ", 's');
%     num_str = num2str(num(i));
%     if num(i) < 10
%         num_str = strcat("000", num_str);
%     end
    fID = fopen(strcat(airfoil_name(i),".pol"),'r');  %loads the file
    A = textscan(fID,'%f %f %f %f %f %f %f', 'HeaderLines', 12);  %skips headers
    fclose(fID);
    alpha{i} = A{:,1};                     %range of alpha
    cl{i} = A{:,2};                        %Coefficient of Lift
    cd{i} = A{:,3};                        %Coefficient of Drag
    cm{i} = A{:,5};                        %Pitching Moment Coefficient
    clcdEfficiency{i} = cl{i}./cd{i};            %drag efficiency
    cont=input("Continue? Y/N ","s");
    if (cont=='Y' || cont=='y')
        i=i+1;
    end
    fclose('all');
end

%% Plotting
figure
hold on
for n=1:i
    plot(alpha{:,n}, cl{:,n})
end
xlabel('alpha');ylabel('CL');title('CL vs. alpha');
legend(airfoil_name, 'Location','southeast')

figure;
hold on
for n=1:i
    plot(alpha{:,n}, cd{:,n})
end
xlabel('alpha');ylabel('CD');title('CD vs. alpha');
legend(airfoil_name)

figure;
hold on
for n=1:i
    plot(alpha{:,n}, cm{:,n})
end
xlabel('alpha');ylabel('CM');title('CM vs. alpha');
legend(airfoil_name)

figure
hold on
for n=1:i
    plot(alpha{:,n}, clcdEfficiency{n})
end
xlabel('alpha');ylabel('CL/CD');title('CL/CD vs. alpha');
legend(airfoil_name)

figure
hold on
for n=1:i
    plot(cl{:,n}, clcdEfficiency{n})
end
xlabel('CL');ylabel('CL/CD');title('CL/CD vs. CL');
legend(airfoil_name,'Location','southwest');
    
    
