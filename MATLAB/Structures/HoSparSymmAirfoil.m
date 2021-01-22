function [Yt] = HoSparSymmAirfoil(x,a4,t,c)
%For the Triton UAS Team
% by: Ryan Beneduce
%   For the determination of max height of spar at desired location
% Source equations found at: airfoiltools.com/airfoil/naca4digit
% Symmetrical airfoils only
% for non-symmetrical airfoils see other function
%inputs of function include: 
%'x'- location along chord length (range from 0-c)
%'a4'- value for closed tail edge= -.1036, normal tail edge= -.1015
%'t'- last two digits of NACA airfoil
%'c'-length of chord
T=t/100;
%Yt calculation
ao=.2969;
a1=-.126;
a2=-.3516;
a3=.2843;
Yt= T*5*c*(ao*(x/c)^.5+a1*x/c+a2*(x/c)^2+a3*(x/c)^3+a4*(x/c)^4);
fprintf('At position x= %d the max spar height is = %d',x,2*Yt)
end

