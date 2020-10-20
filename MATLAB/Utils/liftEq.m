%%Author: Andrew Fletcher
%%Email: afletche@ucsd.edu
%%This program is a MATLAB version of L=W=1/2*rho*velocity^2 *S * Cl
%%Note for emphasis, lift and weight can be interchanged for level flight

%@param hasLift boolean to see if a lift was input
%@param lift the lift that the wings are generating
%@param hasRho boolean to see if a rho was input
%@param rho the density of the fluid (probably air, 1.225 at sea level)
%@param hasVel boolean to see if a velocity was input
%@param velocity the velocity of the airflow (or plane in our case)
%@param hasS boolean to see if an S was input
%@param s wing area
%@param hasCl boolean to see if a Cl was input
%@param cl the lift coefficient
%@return lift see param
%@return rho see param
%@return velocity see param
%@return s see param
%@return cl see param
function [lift, rho, velocity, s, cl] = liftEq(hasLift,lift, hasRho,rho,...
    hasVel,velocity, hasS,s, hasCl,cl)

sum = hasLift + hasRho + hasVel + hasS + hasCl;
% if sum < 4
%     disp('you must already have all the variables but one');
%     disp('One equation, one unknown');
if ~hasLift
    lift = (1/2*rho*velocity^2)*s*cl;   %q = 1/2*rho*velocity^2
elseif ~hasRho
    rho = lift*2/velocity^2/s/cl;
elseif ~hasVel
    velocity = sqrt(lift/s/cl*2/rho);
elseif ~hasS
    s = lift/(1/2*rho*velocity^2)/cl;
elseif ~hasCl
    cl = lift/(1/2*rho*velocity^2)/s;
end

