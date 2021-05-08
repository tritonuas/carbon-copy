A = 0.02; 
B = 0.2; 
C = 0.15; 
D = .25; 
dX = 0.01;
inc = 0.5;
modairfoil_filename = 'C:\\TUAS\\carbon-copy\\MATLAB\\Airfoil\\modairfoil.dat';

Airfoil_Specs = [A B C D;...         %a
                 A+dX B C D;...      %b
                 A B+dX C D;...      %c
                 A B C+dX D;...      %d
                 A B C D+dX;...      %e
                 A+2*dX B C D;...    %f
                 A+dX B+dX C D;...   %g
                 A+dX B C+dX D;...   %h
                 A+dX B C D+dX;...   %i
                 A B+2*dX C D;...    %j
                 A B+dX C+dX D;...   %k
                 A B+dX C D+dX;...   %l
                 A B C+2*dX D;...    %m
                 A B C+dX D+dX;...   %n
                 A B C D+2*dX];      %o
             
for i = 1:15 %For each airfoil: find score
    disp(Airfoil_Specs(i,:));
            modify_airfoil(Airfoil_Specs(i,:),modairfoil_filename); % Create coordinates of airfoil

%             [ClCdScore(i), alpha] = GetClCdScore(Cl,modairfoil_filename); % Find ClCd Score and starting alpha

            [StallScore(i), alphas, Aclcds] = GetStallScore(alpha,modairfoil_filename); % Find stall score 
            [BFStallScore(i), maxAlpha, maxClCd] = BFGetStallScore(modairfoil_filename,inc);
end