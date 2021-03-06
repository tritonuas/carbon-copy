%%This is the main script (executable) for the plane design program
%%As of now, this takes in inputs of wingspan, weight, and initial guesses
%%for the sizing of the parts, and outputs the ideal lift coefficient 
%%and wing reference area to maximize cl/cd and therefore maximize range 
%%and minimize drag for cruise conditions.
%%There are currently 2 approaches to this. One is an analyic method that
%%is based off of the method described in Anderson's 155A lecture slides of
%%setting cdi=cd0. This has been modified such that instead of solving for
%%velocity, it takes in a range of velocity as iterative inputs, and
%%solves for the ideal wing area for each velocity.
%%is done using the minDragEq method. The other method is a purely 
%%iterative method for the purpose of checking.
%%Both processes iterate over velocity to try to find an ideal velocity
%%for the weight. I am currently trying to use the analytic method for the
%%solution.
%%Author: Andrew Fletcher
%%Email: afletche@ucsd.edu

%%Sources: See functions for specific sources. 

%%Concerns: a very high sensitivity to changes in cd0, and
%%the input velocity. This is a concern because variance in manufacturing
%%and variances in flight are likely to be great enough to drastically vary
%%the optimal design according to this program.
%%Also, there are discontinuities arising from the transition point for
%%some reason. This can be seen in the iterative solution.

%%Note: If you see a function being called that you aren't familiar with, I
%%probably made the function. If the function is trivial, then it should be
%%saved under the Utils folder. If not, then it should be in the main
%%folder.

%% General
clear all;
clc;

format compact
format shortg;
addpath("Utils")
addpath("Weight");
addpath("Structures");
addpath("tail");
addpath("Airfoil");
%% Constants
g = 9.81;                   %sea level Earth gravity
density = 1.225;            %sea level air density
viscosity = 1.789*10^-5;    %sea level air viscosity

%% Configuration
%subsonic plane
sweepAngle = 0;     %In the future, I can make this a function of velocity
taperRatio = getTaper(sweepAngle);  %see function for source of eq
                                    %taper should equal 0.45 for subsonic
%% Inputs
% I might make this whole program a function in the future
weight = 135;
lift = weight;      %want to design this for cruise conditions
wingSpan = 3.65;

velocityMin = 5;   %This is a range because both solutions iterate over
velocityMax = 30;
numSteps = 250;
velocity = linspace(velocityMin, velocityMax, numSteps);

sMin = 0.1;     %the iterative solution iterates over these S.
sMax = 2;     %this does not affect the analytic solution
numSteps = 250;
S = linspace(sMin, sMax, numSteps);
SMaster = S;
AR = wingSpan^2./S;  %Getting ARs for induced drag calculations
ARMaster = AR;

%Fuselage Requirements
saFuse = 1;      %1 for carbon copy
lenFuse = 1.2;  %1.2 for carbon copy

%Tail Boom Requirement
tail_boom_radius = 0.0762/2;    %3 inch diameter converted to radius in m
%Note: This is only cause that's the size of the pole we have. This can be
%made a design variable in the future

%%Preparing for the getZeroLiftDrag function by stating that we do not 
%%have these variables (guesses will be given from within the function)
%%-1 means we do not have the information

saNose = -1;
lenNose = -1;
tail_area_h = -1;
cHS = -1;
tail_area_v = -1;
cVS = -1;
cTail = -1;
saTailBoom = -1;
lenTailBoom = -1;

%% Desired Parameters/Constraints

%%CONSTRAINT: MAX CL CRUISE
maxClCruise = 1.8;  %see iterative solution for implementation
maxLoadFactorStall = 1.8;     %input max load factor for takeoff/landing
maxLoadFactorTurns = 1.667; %input max load factor for turns

%%CONSTRIANT: ROOT CHORD LENGTH (STRUCTURAL)
%This is for structural reasons
rootChordReq = 0.3048;  %carbon copy is 0.3048 this is 12 inches
  %That value is based on intuition/history and should be iterated on
tipChordReq = rootChordReq*taperRatio;
avgChordReq = (rootChordReq+tipChordReq)/2;
sMinStruct = avgChordReq*wingSpan;
if sMinStruct > sMin
    sMin = sMinStruct;
    S = linspace(sMin, sMax, numSteps);
    SMaster = S;
    AR = wingSpan^2./S;  %Getting ARs for induced drag calculations
    ARMaster = AR;
end

%%CONSTRAINT: STALL SPEED
stallSpeed = 30;
% L = 1/2*density*stallSpeed^2*S*clMax    nMax = clMax/cl    
%stallSpeed = velocity/sqrt(nMax)   nMax = (velocity/stallSpeed)^2
nStall = (velocityMax/stallSpeed)^2;
if nStall > maxLoadFactorStall
   nStall = maxLoadFactorStall;
   velocityMax = stallSpeed*sqrt(nStall);
   velocity = linspace(velocityMin, velocityMax, numSteps);
end

%%CONSTRAINT: TURNING RADIUS
radius = 99;    %want this radius max 33.98m
%radius = v^2/(g*sqrt(n^2-1))
n = sqrt((velocityMax^2/radius/g)^2 + 1);%load factor req for desired turns
if n > maxLoadFactorTurns
    n = maxLoadFactorTurns;
    velocityMax = sqrt(radius*g*sqrt(n^2-1));  %velocity req for turns
    velocity = linspace(velocityMin, velocityMax, numSteps);
end


%%Add climb rate constraint TODO
%%PROBLEM: This requires a P/W, but it's hard to have information on thrust
climbAngleReq = 10;     %in degrees
G = sind(climbAngleReq);
%Continue to work on this

%%CONSTRAINT: MISSION TIME
%Mission should be completed within 20 minutes.
timeLimit = 9999;
SFTime = 1.05;   %safety factor
timeLimit = timeLimit/SFTime;
timeLimit = timeLimit*60;   %time limit in seconds
%Consult with Garrett to determine whether we should take the whole time or
%cut out time for making second passes.
totalTravelDist = 17000;
additionalWaypointDist = 10000;
totalTravelDist = totalTravelDist + additionalWaypointDist;
velocityReq = totalTravelDist/(timeLimit);
if velocityReq > velocityMin
    velocityMin = velocityReq;
    velocity = linspace(velocityMin, velocityMax, numSteps);
end

%%CONTRAINT: TAIL SIZING
%This is based off of Anderson's MAE 155A slides(see stability, last slide)
C_HT = 0.8;       %Horizontal tail volume coefficient
C_VT = 0.04;     %Vertical tail volume coefficient


%% Function definition
% outputs: function [maxClOverCd_sol, vel_sol, best_cl_sol, best_S_sol, wing_loading, best_chord_sol, best_AR_sol, nIter, stallSpeedIter] ...
% inputs:  = subfunction_name(AR, ARMaster, G, S, SFTime, SMaster, additionalWaypointDist, avgChordReq, cTail, climbAngleReq, density, g, lenFuse, lenNose, lenTailBoom, lift, maxClCruise, maxLoadFactorStall, maxLoadFactorTurns, n, nStall, numSteps, radius, rootChordReq, saFuse, sMax, sMin, sMinStruct, saNose, sTail, saTailBoom, stallSpeed, sweepAngle, taperRatio, timeLimit, tipChordReq, totalTravelDist, velocity, velocityMax, velocityMin, velocityReq, viscosity, weight, wingSpan) 

%Pick Method
methods = ["DBI" "Step_Iterative"];
sol = 0; %has a solution been found?)

if find(methods == "DBI")
    % run the subfunction dbi_sol [with DBI for velocity v and wing area S]
    [maxClOverCd_sol, vel_sol, best_cl_sol, best_S_sol, wing_loading, best_chord_sol, best_AR_sol, nIter, stallSpeedIter] ...
    = dbi_sol(density, g, lenFuse, lenNose, lenTailBoom, lift, maxClCruise, maxLoadFactorStall, maxLoadFactorTurns, radius, saFuse, sMinStruct, sMax, saNose, tail_area_h, cHS, tail_area_v, cVS, saTailBoom, tail_boom_radius, sweepAngle, taperRatio, viscosity, wingSpan, C_HT, C_VT);
    sol = 1;
end

if find(methods == "Step_Iterative")
    % run the subfunction iter_sol [same as before]
    [maxClOverCd_sol, vel_sol, best_cl_sol, best_S_sol, wing_loading, best_chord_sol, best_AR_sol, nIter, stallSpeedIter] ...
    = step_iter_sol(ARMaster, S, SMaster, density, g, lenFuse, lenNose, lenTailBoom, lift, maxClCruise, maxLoadFactorStall, maxLoadFactorTurns, radius, saFuse, saNose, tail_area_h, cHS, tail_area_v, cVS, saTailBoom, sweepAngle, taperRatio, velocity, viscosity, wingSpan, C_HT, C_VT);
    sol = 1;
end
if sol == 0
    disp("Method Selection Invalid. Valid options are: 'DBI', 'Iterative', and 'Analytical'")
end

%% Derivative based iteration
function [clOverCd, v, cl, S, wing_loading, chord, AR, load_factor, stall_speed] ...
    = dbi_sol(density, g, lenFuse, lenNose, lenTailBoom, lift, maxClCruise, maxLoadFactorStall, maxLoadFactorTurns, radius, saFuse, sMinStruct, sMax, saNose, tail_area_h, cHS, tail_area_v, cVS, saTailBoom, tail_boom_radius, sweepAngle, taperRatio, viscosity, wingSpan, C_HT, C_VT) 
%TODO still need to finish validating and cleaning up stuff that is not
%needed, and finish making the output print statements
%TODO think about this from a controller standpoint, and how to add Kd and
%Ki terms if I need to.
%declaring that I have these variables for use of the lift equation
hasLift = 1;
hasDensity = 1;
hasVel = 1;

hasS = 1;       %declaring that I have S and want to calculate for Cl
hasCl = 0;      %This is in preperation for the lift equation

S = 5;
Kp_S = 1;
iterNum = 0;
derivative_cl_over_cd = 1;   %to run the loop
delta_S = 0.01; %first step
while abs(derivative_cl_over_cd) > 0.001
    iterNum = iterNum + 1;
%     
%     while abs(derivative_S) > 0.001
%         iterNum = iterNum + 1;

    %Placing a minimum S for structural 
    if S(iterNum) < sMinStruct
         S(iterNum) = sMinStruct;
    end 
    %Placing a maximum S  
    if S(iterNum) > sMax
         S(iterNum) = sMax;
    end
    
    % calculate parameters based on S (chord, AR)
    chord = S(iterNum)/wingSpan;   %from max cl checker
    AR = wingSpan^2/S(iterNum); 

    %calculate weight
    weight = getWeight(S(iterNum));       
    lift = weight;

    v_guess = 20;
    v(iterNum) = 10;
    while abs(v_guess - v(iterNum)) > 0.1
        v_guess = v(iterNum);
        
        saWing = S(iterNum)*2;
%         saTail = sTail*2;
        saHS = tail_area_h*2;
        saVS = tail_area_v*2;
        %get the zero-lift drag coeff for this S and velocity
        cd0(iterNum) = getZeroLiftDrag(density, viscosity, v(iterNum), ...
                   S(iterNum), saWing, chord, saFuse,lenFuse, saNose,lenNose,...
                   saHS,cHS, saVS,cVS, saTailBoom,lenTailBoom);

        %calculate velocity and lift coefficient
        v(iterNum) = -1;
        [S(iterNum), v(iterNum), cl(iterNum)] = minDragEq(v(iterNum),...
          density, S(iterNum), cd0(iterNum), wingSpan, lift, taperRatio, sweepAngle);
        saWing = S(iterNum)*2;
      
        % calculate tail areas and tail boom length (and therefore area)
        [tail_area_h, tail_area_v, lenTailBoom] = find_tail_size(wingSpan, ...
        S(end), saWing, chord(end), C_HT, C_VT, density, viscosity, v(end));
        sTail = tail_area_h + tail_area_v;
        saTailBoom = lenTailBoom*2*pi*tail_boom_radius;
    end
  
    q = 1/2*density*v(iterNum)^2;
    
    %Placing a maximum lift coefficient due to airfoil constraints
    if cl(iterNum) > maxClCruise
         cl(iterNum) = maxClCruise;
         hasS = 0;
         hasCl = 1;
         %Overriding the S to meet the new constraint with a fixed cl
         [lift,density,v(iterNum),S(iterNum), cl(iterNum)] = liftEq(hasLift,lift,...
         hasDensity,density,hasVel,v(iterNum),hasS,S(iterNum),hasCl,cl(iterNum));
         %Resetting design variables for future iterations
         hasS = 1;
         hasCl = 0;
         chord(iterNum) = S(iterNum)/wingSpan;   %updating chord for new S
         AR = wingSpan^2/S(iterNum);    %updating asepct ratio for new S
    end

    %get induced drag so we can later get cd
    k = getK(AR, taperRatio, sweepAngle);
    cdi(iterNum) = cl(iterNum)^2*k;
    %get cd so we can later get cl/cd
    cd(iterNum) = cd0(iterNum) + cdi(iterNum);
    %get cl/cd for performance comparison
    clOverCd(iterNum) = cl(iterNum)/cd(iterNum);

    if iterNum == 1
        %If initial S > sMax (and could not check smaller values)
        if S(iterNum) == sMax
            delta_S = -delta_S;
        end
        S(iterNum+1) = S(iterNum) + delta_S;
    else
        d_cl_over_cd = clOverCd(iterNum) - clOverCd(iterNum-1);
        derivative_cl_over_cd(iterNum) = d_cl_over_cd/delta_S;
        delta_S = Kp_S*derivative_cl_over_cd(iterNum);
        S(iterNum+1) = S(iterNum) + delta_S;
    end
end
S = S(1:end-1);

figure(1);
plot(2:iterNum, clOverCd(2:end), '-bo');
xlabel("Iteration Number"); ylabel("CL/CD"); 
title("Derivative-Based Iteration: Cl/CD vs. Iteration Number");
figure(2);
plot(2:iterNum, S(2:end), '-bo');
xlabel("Iteration Number"); ylabel("Wing Area (m^2)"); 
title("Derivative-Based Iteration: Wing Area vs. Iteration Number");

%update load factor
load_factor = sqrt((v(iterNum)^2/radius/g)^2 + 1);
%update stall speed
stall_speed = v(end)/sqrt(maxLoadFactorStall);

wing_loading = lift/S(end);

disp("----------------------------");
disp("DBI_ANA SOLUTION");
disp("----------------------------");
disp("Max cl/cd: " + clOverCd(end));
disp("Velocity: " + v(end));
disp("cl: " + cl(end));
disp("S: " + S(end));
disp("Wing loading: " + wing_loading);
disp("Chord: " + chord);
disp("Aspect ratio: " + AR);
disp("Load factor on turns: " + load_factor);
disp("Horizontal Tail Area: " + tail_area_h); 
disp("Vertical Tail Area: " + tail_area_v); 
disp("Tail Boom Length: " + lenTailBoom); 
if load_factor > maxLoadFactorTurns + 0.001
    disp("Invalid solution! Load factor constraint not satisfied.");
end
disp("Stall speed: " + stall_speed);
% 
end
% 

%% Stepwise,Nested Iterative solution
function [maxClOverCd_sol, vel_sol, best_cl_sol, best_S_sol, wing_loading, best_chord_sol, best_AR_sol, nIter, stallSpeedIter] ...
    = step_iter_sol(ARMaster, S, SMaster, density, g, lenFuse, lenNose, lenTailBoom, lift, maxClCruise, maxLoadFactorStall, maxLoadFactorTurns, radius, saFuse, saNose, tail_area_h, cHS, tail_area_v, cVS, saTailBoom, sweepAngle, taperRatio, velocity, viscosity, wingSpan, C_HT, C_VT) 
%declaring that I have these variables for use of the lift equation
hasLift = 1;
hasDensity = 1;
hasVel = 1;

hasS = 1;       %declaring that I have S and want to calculate for Cl
hasCl = 0;      %This is in preperation for the lift equation

%%preallocation for inner loop (S)
cl = S;
cd = S;
cd0 = S;
cdi = S;
clOverCd = S;
drag = S;
liftGen = S;

%preallocation for outer loop (velocity)
maxClOverCdIterate = velocity;
bestS = velocity;
bestchord = velocity;
bestAR = velocity;
bestcd0 = velocity;
bestcdi = velocity;
bestcd = velocity;
bestcl = velocity;
bestDrag = velocity;
bestLift = velocity;
zeroLiftDrag = velocity;
inducedDrag = velocity;

maxOfMaxClOverCd = 0;
for j = 1:length(velocity)  %for every velocity...
    maxClOverCd = 0;
    S = SMaster;  %This is to remove any changes from proir velocities
    chord = SMaster/wingSpan;
    AR = ARMaster;
    
    for i = 1:length(S) %for every wing area...
        weight = getWeight(S(i));       
        lift = weight; 
        
        %get the lift coefficient necessary for this S and velocity
        [lift, density, velocity(j), S(i), cl(i)] = liftEq(hasLift,lift,...
           hasDensity,density, hasVel,velocity(j), hasS,S(i), hasCl,cl(i));
        %Placing a maximum lift coefficient due to airfoil constraints
        if cl(i) > maxClCruise
            cl(i) = maxClCruise;
            hasS = 0;
            hasCl = 1;
            %Overriding the S to meet the new constraint with a fixed cl
            [lift,density,velocity(j),S(i),cl(i)] = liftEq(hasLift,lift,...
            hasDensity,density,hasVel,velocity(j),hasS,S(i),hasCl,cl(i));
            %Resetting design variables for future iterations
            hasS = 1;
            hasCl = 0;
            chord(i) = S(i)/wingSpan;   %updating chord for new S
            AR(i) = wingSpan^2/S(i);    %updating asepct ratio for new S
        end

        saWing = S(i)*2;
        saHS = tail_area_h*2;
        saVS = tail_area_v*2;
        cd0(i) = getZeroLiftDrag(density, viscosity, velocity(j), ...
                   S(i), saWing, chord(i), saFuse,lenFuse, saNose,lenNose,...
                   saHS, cHS, saVS, cVS, saTailBoom,lenTailBoom);
               
        %get induced drag so we can later get cd
        k = getK(AR(i), taperRatio, sweepAngle);
        cdi(i) = cl(i)^2*k;
        %get cd so we can later get cl/cd
        cd(i) = cd0(i) + cdi(i);
        %get cl/cd for performance comparison
        clOverCd(i) = cl(i)/cd(i);
        %get max cl/cd for best S for each velocity
        if clOverCd(i) >= maxClOverCd
            maxClOverCd = clOverCd(i);
            maxIndex = i;
        end
        %Get drag for plotting
        drag(i) = dragEq(0,0, hasDensity,density, hasVel,velocity(j),...
            hasS,S(i), 1,cd(i));
        
    end
    
    %get characteristics of each local best performance
    maxClOverCdIterate(j) = maxClOverCd;
    bestS(j) = S(maxIndex);
    bestchord(j) = bestS(j)/wingSpan;
    bestAR(j) = wingSpan^2/bestS(j);
    bestcd0(j) = cd0(maxIndex);
    zeroLiftDrag(j) = 1/2*density*velocity(j)^2*bestS(j)*bestcd0(j); %for plotting purposes
    bestcdi(j) = cdi(maxIndex);
    inducedDrag(j) = 1/2*density*velocity(j)^2*bestS(j)*bestcdi(j); %for plotting purposes
    bestcd(j) = cd(maxIndex);
    bestcl(j) = cl(maxIndex);
    bestDrag(j) = drag(maxIndex);
    
    %get best performance for a sea level airframe of this weight
    if maxClOverCdIterate(j) > maxOfMaxClOverCd
        maxOfMaxClOverCd = maxClOverCdIterate(j);
        maxOfMaxIndex = j;  %can use to get characteristics of best
    end
end

%update load factor
nIter = sqrt((velocity(maxOfMaxIndex)^2/radius/g)^2 + 1);
%Update stall speed
stallSpeedIter = velocity(maxOfMaxIndex)/sqrt(maxLoadFactorStall);

%%Plots for visualization
figure(101)
plot(velocity, bestDrag, '-k', ...
    velocity, zeroLiftDrag, '-r', velocity, inducedDrag, '-m')
xlabel("Velocity (m/s)"); ylabel("Drag (N)"); 
title("Iterative: Minimum Drag vs. Velocity");
legend("total drag","zero lift drag","induced drag");
figure(102)
plot(velocity, maxClOverCdIterate, '-b')
xlabel("Velocity (m/s)"); ylabel("CL/CD"); 
title("Iterative: Maximum CL/CD vs. Velocity");

maxClOverCd_sol = maxOfMaxClOverCd;
vel_sol = velocity(maxOfMaxIndex);
best_cl_sol = bestcl(maxOfMaxIndex);
best_S_sol = bestS(maxOfMaxIndex);

%Update weight and lift
weight = getWeight(best_S_sol);
lift = weight;

wing_loading = lift/bestS(maxOfMaxIndex);
best_chord_sol = bestchord(maxOfMaxIndex);
best_AR_sol = bestAR(maxOfMaxIndex);

[tail_area_h, tail_area_v, tail_boom_length] = find_tail_size(wingSpan, ...
    best_S_sol, best_S_sol*2, best_chord_sol, C_HT, C_VT, density, viscosity, vel_sol);

disp("----------------------------");
disp("ITERATIVE SOLUTION");
disp("----------------------------");
disp("Max cl/cd: " + maxClOverCd_sol);
disp("Velocity: " + vel_sol);
disp("cl: " + best_cl_sol);
disp("S: " + best_S_sol);
disp("Wing loading: " + wing_loading);
disp("Chord: " + best_chord_sol);
disp("Aspect ratio: " + best_AR_sol);
disp("Load factor on turns: " + nIter);
disp("Horizontal Tail Area: " + tail_area_h); 
disp("Vertical Tail Area: " + tail_area_v); 
disp("Tail Boom Length: " + tail_boom_length); 
if nIter > maxLoadFactorTurns + 0.001
    disp("Invalid solution! Load factor constraint not satisfied.");
end
disp("Stall speed: " + stallSpeedIter);
end

