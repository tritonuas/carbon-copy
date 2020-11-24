%%This function gets the T transformation matrix to convert the
%%from the X,Y coordinate system to the 11,22,12 (planar case)
%%Author: Andrew Fletcher
%%Email: afletche@ucsd.edu
%@param theta the angle of the lamina
%@param stressOrStrain a string of 'stress' or 'strain' to denote which one
%@param radOrDeg a string of 'rad' or 'deg' tp denote units of theta
%@return T_sigma the stress transformation matrix
function T = get_T_planar(theta, stressOrStrain, radOrDeg)

%defining the units of theta
if radOrDeg == "rad"
    c = cos(theta);
    s = sin(theta);
elseif radOrDeg == "deg"
    c = cosd(theta);
    s = sind(theta);
else
    disp("Please input 'rad' for radians or 'deg' for degrees");
end

%defining whether we are solving for T_sigma or T_epsilon
if stressOrStrain == "stress"
    T = [c^2  s^2  2*s*c;    ...
         s^2  c^2 -2*s*c;    ...
         -s*c s*c (c^2-s^2)];
elseif stressOrStrain == "strain"
    T = [c^2    s^2    s*c;    ...
         s^2    c^2   -s*c;    ...
         -2*s*c 2*s*c (c^2-s^2)];
else
    disp("Please input 'stress' for T_sigma and 'strain' for T_epsilon");
end