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

%%Author: Andrew Fletcher
%%Email: afletche@ucsd.edu
%%Note for emphasis, lift and weight can be interchanged for level flight

%%Source for equation: Professor Anderson (MAE 155A lecture slides)
    %Specifically the slides called "Level" for level flight

%%Note: This equation is originally configured to solve for the velocity,
%when given an S, and everything else, but I am planning on using it to
%solve for S given velocity.

%%Assumptions: L=W, drag is minimized when Cl/Cd is maximized (true
%statement), and Cl/Cd is maximized when (del/del(Cl))(Cl/Cd) = 0; This
%happens when cl is taken as sqrt(cd0/k)

%%Concerns: only the same concerns as from the analytic solution. Main
%concern in here comes from the uncertainty in the getOswaldEff function.

function [s, velocity, cl, k] = minDragEq(velocity, density, s, cd0, ...
    wingSpan, lift, taperRatio, sweepAngle)

if density == -1
    disp("Default density is 1.225 kg/m^3 (sea level air)");
end
if taperRatio == -1
    taperRatio = 0.45;
    disp("Default taper ratio is 0.45 (good for subsonic)");
end
if sweepAngle == -1
    sweepAngle = 0;
    disp("Default taper ratio is 0 degrees sweep (good for subsonic)");
end
if wingSpan == -1
    wingSpan = 3;
    disp("Default wing span is 3m. Please input a wing span");
end
if cd0 == -1
    cd0 = 0.05;
    disp("Default cd0 is 0.05 which is totally arbitrary." +...
        "Please input a cd0");
end
if lift == -1
    lift = 135;
    disp("Default lift=weight is 135 N. Please input a weight");
end


if s == -1 && velocity == -1
    disp("Please input either a velocity or a wing area.")
elseif s == -1
    eGuess = 0.9;
    eReal = eGuess;
    eGuess = 0; %so that it runs the loop
    counter = 0;
    while(abs(eReal - eGuess) > 10^-7)
        eGuess = eReal;
        s = 4*lift^2/(velocity^4*density^2*wingSpan^2*pi*eGuess*cd0);
        ar = wingSpan^2/s;
        eReal = getOswaldEff(taperRatio, ar, sweepAngle);
        if counter > 1000
            disp("stuck in oscillation between eGuess and eReal in minDragEq");
            break
        end
    end
    k = 1/(pi*eReal*ar);
    cl = sqrt(cd0/k);   %the condition that builds this equation
elseif velocity == -1
    ar = wingSpan^2/s;
    k = getK(ar, taperRatio, sweepAngle);
    velocity = sqrt(2*lift/(density*s)*sqrt(k/cd0));
    cl = sqrt(cd0/k);   %the condition that builds this equation
end
