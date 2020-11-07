%% XFOIL 2 Matlab
% Given a folder of airfoil points, uses XFOIL to sims on the airfoil
% through matlab, saves data in a pol file
% Kaelan Tan 11/2020

%%
% Set up the file names
A=dir('airfoil_database');  % include xfoil in same folder as airfoil dat files
names={A.name}; names(1)=[]; names(1)=[];  % First two files arent airfoils

% Code we want to run on Xfoil
for i=1:length(names)
    fname=char(names(i));   % formatting for ease
    name=fname(1:end-4);
    named(i)=string(name);
    fid = fopen('xfoil_input.txt','w');   % create the inputs
    fprintf(fid,['load '+named(i)+'.dat\n']);   % commands
    fprintf(fid,[named(i)+'.dat\n']);   % goes through the names of the files in Airoils folder
%     fprintf(fid,['load '+named(i)+'\n']);   % commands
%     fprintf(fid,[named(i)+'\n']);   % goes through the names of the files in Airoils folder
    fprintf(fid,'pane\n');   % makes the airfoils nice
    fprintf(fid,['oper\n']);fprintf(fid,['visc\n']);   % makes visc analysis w/ Re=4E5
    fprintf(fid,['4e5\n']);fprintf(fid,['pacc\n']);
    fprintf(fid,[named(i)+'.pol\n\n']);
    fprintf(fid,'iter\n');fprintf(fid,'25\n'); % aint nobody got time fo 50
    fprintf(fid,'aseq 0 20 0.5\n');% alpha 0 to 20 deg in 0.5 deg increment
    fprintf(fid,'pacc\n');

    cmd = 'xfoil.exe < xfoil_input.txt';   % running on xfoil
    [status,result] = system(cmd);

    fclose('all');
%     delete('xfoil_input.txt');
end
