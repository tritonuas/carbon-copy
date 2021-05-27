A = 0.01; 
B = 0.2; 
C = 0.15; 
D = .25; 
Cl = 0.2;
dX = 0.01;
alpha = 0;
modairfoil_filename = 'C:\\TUAS\\carbon-copy\\MATLAB\\Airfoil\\modairfoil.dat';

Airfoil_Specs = [A B C D;...         %a
                 A+dX B C D;...      %b
                 A B+dX C D;...      %c             %two point error in GetStallScore
                 A B C+dX D;...      %d
                 A B C D+dX;...      %e
                 A+2*dX B C D;...    %f             %two point error in GetStallScore
                 A+dX B+dX C D;...   %g
                 A+dX B C+dX D;...   %h
                 A+dX B C D+dX;...   %i
                 A B+2*dX C D;...    %j
                 A B+dX C+dX D;...   %k             %Converges to local max, not absolute max
                 A B+dX C D+dX;...   %l
                 A B C+2*dX D;...    %m             %Good example of converged yet bouncing
                 A B C+dX D+dX;...   %n
                 A B C D+2*dX];      %o             %bounces between two alphas
             
               Airfoil_Specs = [A B+dX C D];
             
for i = 1:height(Airfoil_Specs) %For each airfoil: find score
    disp(Airfoil_Specs(i,:));
            modify_airfoil(Airfoil_Specs(i,:),modairfoil_filename); % Create coordinates of airfoil
            [ClCdScore(i), alpha] = GetClCdScore(Cl,modairfoil_filename);
%             disp(alpha);
            [StallScore(i), alphas, Aclcds] = GetStallScore(alpha,modairfoil_filename); % Find stall score 
%             [BFStallScore(i), maxAlpha, maxClCd] = BFGetStallScore(modairfoil_filename);
end
disp(StallScore);