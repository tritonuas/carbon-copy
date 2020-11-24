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
%addpath("Airfoil");
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
sFuse = 1;      %1 for carbon copy
lenFuse = 1.2;  %1.2 for carbon copy

%%Preparing for the getZeroLiftDrag function by stating that we do not 
%%have these variables (guesses will be given from within the function)
%%-1 means we do not have the information
sNose = -1;
lenNose = -1;
sTail = -1;
cTail = -1;
sTailBoom = -1;
lenTailBoom = -1;

%% Desired Parameters

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
radius = 99;    %want this radius
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
C_HT = 1;       %Horizontal tail volume coefficient
C_VT = 0.1;     %Vertical tail volume coefficient


%% Function definition
% outputs: function [maxClOverCd_sol, vel_sol, best_cl_sol, best_S_sol, wing_loading, best_chord_sol, best_AR_sol, nIter, stallSpeedIter] ...
% inputs:  = subfunction_name(AR, ARMaster, G, S, SFTime, SMaster, additionalWaypointDist, avgChordReq, cTail, climbAngleReq, density, g, lenFuse, lenNose, lenTailBoom, lift, maxClCruise, maxLoadFactorStall, maxLoadFactorTurns, n, nStall, numSteps, radius, rootChordReq, sFuse, sMax, sMin, sMinStruct, sNose, sTail, sTailBoom, stallSpeed, sweepAngle, taperRatio, timeLimit, tipChordReq, totalTravelDist, velocity, velocityMax, velocityMin, velocityReq, viscosity, weight, wingSpan) 

%Pick Method
methods = ["DBI" "DBI_it" "Iterative" "Analytical"];
sol = 0; %has a solution been found?)

if find(methods == "DBI")
    % run the subfunction dbi_sol [with DBI for velocity v and wing area S]
    [maxClOverCd_sol, vel_sol, best_cl_sol, best_S_sol, wing_loading, best_chord_sol, best_AR_sol, nIter, stallSpeedIter] ...
    = dbi_sol(cTail, density, g, lenFuse, lenNose, lenTailBoom, lift, maxClCruise, maxLoadFactorStall, maxLoadFactorTurns, radius, sFuse, sMinStruct, sMax, sNose, sTail, sTailBoom, sweepAngle, taperRatio, viscosity, wingSpan, C_HT, C_VT);
    sol = 1;
end

if find(methods == "DBI_it")
    % run the subfunction dbi_sol [same as before, iterative wing area]
    [maxClOverCd_sol, vel_sol, best_cl_sol, best_S_sol, wing_loading, best_chord_sol, best_AR_sol, nIter, stallSpeedIter] ...
    = dbi_itS_sol(ARMaster, S, SMaster, cTail, density, g, lenFuse, lenNose, lenTailBoom, lift, maxClCruise, maxLoadFactorStall, maxLoadFactorTurns, radius, sFuse, sNose, sTail, sTailBoom, sweepAngle, taperRatio, viscosity, wingSpan, C_HT, C_VT);
    sol = 1;
end

if find(methods == "Iterative")
    % run the subfunction iter_sol [same as before]
    [maxClOverCd_sol, vel_sol, best_cl_sol, best_S_sol, wing_loading, best_chord_sol, best_AR_sol, nIter, stallSpeedIter] ...
    = iter_sol(ARMaster, S, SMaster, cTail, density, g, lenFuse, lenNose, lenTailBoom, lift, maxClCruise, maxLoadFactorStall, maxLoadFactorTurns, radius, sFuse, sNose, sTail, sTailBoom, sweepAngle, taperRatio, velocity, viscosity, wingSpan, C_HT, C_VT);
    sol = 1;
end
if find(methods == "Analytical")
    % run the outdated subfunction anal_sol
    % no comparison graphs
    [maxClOverCd_sol, vel_sol, best_cl_sol, best_S_sol, wing_loading, best_chord_sol, best_AR_sol, nIter, stallSpeedIter] ...
    = anal_sol(cTail, density, g, lenFuse, lenNose, lenTailBoom, lift, maxLoadFactorStall, maxLoadFactorTurns, radius, sFuse, sNose, sTail, sTailBoom, sweepAngle, taperRatio, velocity, viscosity, wingSpan, C_HT, C_VT);
    sol = 1;
end
if sol == 0
    disp("Method Selection Invalid. Valid options are: 'DBI', 'Iterative', and 'Analytical'")
end

%% Derivative based iteration solution, iterative wing area
function [maxClOverCd_sol, vel_sol, best_cl_sol, best_S_sol, wing_loading, best_chord_sol, best_AR_sol, nIter, stallSpeedIter] ...
    = dbi_sol(cTail, density, g, lenFuse, lenNose, lenTailBoom, lift, maxClCruise, maxLoadFactorStall, maxLoadFactorTurns, radius, sFuse, sMinStruct, sMax, sNose, sTail, sTailBoom, sweepAngle, taperRatio, viscosity, wingSpan, C_HT, C_VT) 
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

zerocol = zeros(250); % random :-\
cl = zerocol(1,:); %%%%%%%%%% Idk why cl needs pre-allocation

v_dbi = 10;
Kp_v = 1.6;
ViterNum = 0;
derivative_v = 1;   %to run the loop
delta_v = 0.01; %first step
while abs(derivative_v) > 0.001
    ViterNum = ViterNum + 1;
    
    % defining for S-iter loop
    S_dbi = 42; % sMinStruct is 0.80658
    Kp_S = 1.6;
    SiterNum = 0;
    derivative_S = 1;   %to run the loop
    delta_S = 0.01; %first step
    
    while abs(derivative_S) > 0.001
        SiterNum = SiterNum + 1;

       %Placing a minimum S for structural 
       if S_dbi(SiterNum) < sMinStruct
            S_dbi(SiterNum) = sMinStruct;
       end 
       
       %Placing a maximum S  
       if S_dbi(SiterNum) > sMax
            S_dbi(SiterNum) = sMax;
       end

       % calculate parameters based on S (chord, AR)
       chord(SiterNum) = S_dbi(SiterNum)/wingSpan;   %from max cl checker
       AR(SiterNum) = wingSpan^2/S_dbi(SiterNum); 
       
       %calculate weight
       weight = getWeight(S_dbi(SiterNum));       
       lift = weight;
       %get the lift coefficient necessary for this S and velocity
       [lift, density, v_dbi(ViterNum), S_dbi(SiterNum), cl(SiterNum)] = liftEq(hasLift,lift,...
          hasDensity,density,hasVel,v_dbi(ViterNum), hasS,S_dbi(SiterNum), hasCl,cl(SiterNum));

       %Placing a maximum lift coefficient due to airfoil constraints
       if cl(SiterNum) > maxClCruise
            cl(SiterNum) = maxClCruise;
            hasS = 0;
            hasCl = 1;
            %Overriding the S to meet the new constraint with a fixed cl
            [lift,density,v_dbi(ViterNum),S_dbi(SiterNum), cl(SiterNum)] = liftEq(hasLift,lift,...
            hasDensity,density,hasVel,v_dbi(ViterNum),hasS,S_dbi(SiterNum),hasCl,cl(SiterNum));
            %Resetting design variables for future iterations
            hasS = 1;
            hasCl = 0;
            chord(SiterNum) = S_dbi(SiterNum)/wingSpan;   %updating chord for new S
            AR(SiterNum) = wingSpan^2/S_dbi(SiterNum);    %updating asepct ratio for new S
       end

        %get the zero-lift drag coeff for this S and velocity
        cd0(SiterNum) = getZeroLiftDrag(density, viscosity, v_dbi(ViterNum), ...
                    S_dbi(SiterNum), chord(SiterNum), sFuse,lenFuse, sNose,lenNose,...
                    sTail,cTail, sTailBoom,lenTailBoom);

        %get induced drag so we can later get cd
        k = getK(AR(SiterNum), taperRatio, sweepAngle);
        cdi(SiterNum) = cl(SiterNum)^2*k;
        %get cd so we can later get cl/cd
        cd(SiterNum) = cd0(SiterNum) + cdi(SiterNum);
        %get cl/cd for performance comparison
        clOverCd(SiterNum) = cl(SiterNum)/cd(SiterNum);
        %get max cl/cd for best S for each velocity
        
        % old maxClOverCd searcher / checker
        maxClOverCd = clOverCd(SiterNum);
        maxIndex = SiterNum;
        
        %Get drag for plotting
        drag(SiterNum) = dragEq(0,0, hasDensity,density, hasVel,v_dbi(ViterNum),...
            hasS,S_dbi(SiterNum), 1,cd(SiterNum));
        
        %get characteristics of each local best performance. in 2d matrix
        %to see iteration?? Could delete
        maxClOverCdIterate_S(SiterNum) = clOverCd(SiterNum); %%%%%%%maxClOverCd;
        bestS(ViterNum, SiterNum) = S_dbi(maxIndex);
        bestchord(ViterNum, SiterNum) = bestS(ViterNum, SiterNum)/wingSpan;
        bestAR(ViterNum, SiterNum) = wingSpan^2/bestS(ViterNum, SiterNum);
        bestcd0(ViterNum, SiterNum) = cd0(maxIndex);
        zeroLiftDrag(ViterNum, SiterNum) = 1/2*density*S_dbi(maxIndex)^2*bestS(ViterNum, SiterNum)*bestcd0(ViterNum, SiterNum); %for plotting purposes
        bestcdi(ViterNum, SiterNum) = cdi(maxIndex);
        inducedDrag(ViterNum, SiterNum) = 1/2*density*S_dbi(maxIndex)^2*bestS(ViterNum, SiterNum)*bestcdi(ViterNum, SiterNum); %for plotting purposes
        bestcd(ViterNum, SiterNum) = cd(maxIndex);
        bestcl(ViterNum, SiterNum) = cl(maxIndex);
        bestDrag(ViterNum, SiterNum) = drag(maxIndex);
        
        if SiterNum == 1
            %If initial S > sMax (and could not check smaller values)
            if S_dbi(SiterNum) == sMax
                S_dbi(SiterNum+1) = S_dbi(SiterNum) - delta_S; % make second guess decrease
            else
                S_dbi(SiterNum+1) = S_dbi(SiterNum) + delta_S;
            end
        else
            delta_cl_over_cd_S = maxClOverCdIterate_S(SiterNum)...
                - maxClOverCdIterate_S(SiterNum-1);           
            derivative_S(SiterNum) = delta_cl_over_cd_S/delta_S;
            delta_S = Kp_S*derivative_S(SiterNum);
            S_dbi(SiterNum+1) = S_dbi(SiterNum) + delta_S;
        end
    end 
    S_dbi = S_dbi(1:end-1); % clip off, not real
    maxClOverCdIterate(ViterNum) = maxClOverCdIterate_S(maxIndex);

    opt_S(ViterNum) = S_dbi(end);
    %get delta_cl_over_cd
    if ViterNum == 1
        v_dbi(ViterNum+1) = v_dbi(ViterNum) + delta_v;
    else
        delta_cl_over_cd = maxClOverCdIterate(ViterNum)...
            - maxClOverCdIterate(ViterNum-1);
        derivative_v(ViterNum) = delta_cl_over_cd/delta_v;
        delta_v = Kp_v*derivative_v(ViterNum);
        v_dbi(ViterNum+1) = v_dbi(ViterNum) + delta_v;
    end
end
v_dbi = v_dbi(1:end-1); %clip off the last value because it's not real

figure(1);
plot(1:ViterNum, maxClOverCdIterate, '-bo');
xlabel("Iteration Number"); ylabel("CL/CD"); 
title("Derivative-Based Iteration: Cl/CD vs. Iteration Number");
figure(2);
plot(1:ViterNum, v_dbi, '-bo');
xlabel("Iteration Number"); ylabel("Velocity (m/s)"); 
title("Derivative-Based Iteration: Velocity vs. Iteration Number");
figure(3); %%%%%%% for your reference, other variables didn't change this dramatically
plot(1:ViterNum, opt_S, '-bo');
xlabel("Iteration Number"); ylabel("Area (m^2)"); 
title("Derivative-Based Iteration: Area vs. Iteration Number");

%update load factor
nIter = sqrt((v_dbi(ViterNum)^2/radius/g)^2 + 1);
%update stall speed
stallSpeedIter = v_dbi(end)/sqrt(maxLoadFactorStall);

maxClOverCd_sol = maxClOverCdIterate(end);
vel_sol = v_dbi(end);
best_cl_sol = bestcl(end, SiterNum);
best_S_sol = bestS(end, SiterNum);
wing_loading = lift/bestS(end, SiterNum);
best_chord_sol = bestchord(end, SiterNum);
best_AR_sol = bestAR(end, SiterNum); %replaced all IterNum with (end) or (end, SiterNum)

% calculate tail areas
[tail_area_h, tail_area_v, tail_boom_length] = find_tail_size(wingSpan, ...
    best_S_sol, best_chord_sol, C_HT, C_VT, density, viscosity, vel_sol);
disp("----------------------------");
disp("DBI SOLUTION");
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
% 
end
% 
%% Derivative based iteration solution (original)
function [maxClOverCd_sol, vel_sol, best_cl_sol, best_S_sol, wing_loading, best_chord_sol, best_AR_sol, nIter, stallSpeedIter] ...
    = dbi_itS_sol(ARMaster, S, SMaster, cTail, density, g, lenFuse, lenNose, lenTailBoom, lift, maxClCruise, maxLoadFactorStall, maxLoadFactorTurns, radius, sFuse, sNose, sTail, sTailBoom, sweepAngle, taperRatio, viscosity, wingSpan, C_HT, C_VT) 
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

%%preallocation for inner loop (S)
cl = S;
cd = S;
cd0 = S;
cdi = S;
clOverCd = S;
drag = S;
liftGen = S;

maxOfMaxClOverCd = 0;
v_dbi = 10;
Kp = 1.6;
iterNum = 0;
derivative = 1;   %to run the loop
delta_v = 0.01; %first step
while abs(derivative) > 0.001
    iterNum = iterNum + 1;
    
    maxClOverCd = 0;
    S = SMaster;  %This is to remove any changes from proir velocities
    chord = SMaster/wingSpan;
    AR = ARMaster;
        
    for i = 1:length(S) %for every wing area...
         weight = getWeight(S(i));       
         lift = weight;
        %get the lift coefficient necessary for this S and velocity
        [lift, density, v_dbi(iterNum), S(i), cl(i)] = liftEq(hasLift,lift,...
           hasDensity,density, hasVel,v_dbi(iterNum), hasS,S(i), hasCl,cl(i));
        %Placing a maximum lift coefficient due to airfoil constraints
        if cl(i) > maxClCruise
            cl(i) = maxClCruise;
            hasS = 0;
            hasCl = 1;
            %Overriding the S to meet the new constraint with a fixed cl
            [lift,density,v_dbi(iterNum),S(i),cl(i)] = liftEq(hasLift,lift,...
            hasDensity,density,hasVel,v_dbi(iterNum),hasS,S(i),hasCl,cl(i));
            %Resetting design variables for future iterations
            hasS = 1;
            hasCl = 0;
            chord(i) = S(i)/wingSpan;   %updating chord for new S
            AR(i) = wingSpan^2/S(i);    %updating asepct ratio for new S
        end

        %get the zero-lift drag coeff for this S and velocity
        cd0(i) = getZeroLiftDrag(density, viscosity, v_dbi(iterNum), ...
                    S(i), chord(i), sFuse,lenFuse, sNose,lenNose,...
                    sTail,cTail, sTailBoom,lenTailBoom);

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
        drag(i) = dragEq(0,0, hasDensity,density, hasVel,v_dbi(iterNum),...
            hasS,S(i), 1,cd(i));
        
    end
    
    %get characteristics of each local best performance
    maxClOverCdIterate(iterNum) = maxClOverCd;
    bestS(iterNum) = S(maxIndex);
    bestchord(iterNum) = bestS(iterNum)/wingSpan;
    bestAR(iterNum) = wingSpan^2/bestS(iterNum);
    bestcd0(iterNum) = cd0(maxIndex);
    zeroLiftDrag(iterNum) = 1/2*density*v_dbi(iterNum)^2*bestS(iterNum)*bestcd0(iterNum); %for plotting purposes
    bestcdi(iterNum) = cdi(maxIndex);
    inducedDrag(iterNum) = 1/2*density*v_dbi(iterNum)^2*bestS(iterNum)*bestcdi(iterNum); %for plotting purposes
    bestcd(iterNum) = cd(maxIndex);
    bestcl(iterNum) = cl(maxIndex);
    bestDrag(iterNum) = drag(maxIndex);
    
    %get delta_cl_over_cd
    if iterNum == 1
        v_dbi(iterNum+1) = v_dbi(iterNum) + delta_v;
    else
        delta_cl_over_cd = maxClOverCdIterate(iterNum)...
            - maxClOverCdIterate(iterNum-1);
        derivative(iterNum) = delta_cl_over_cd/delta_v;
        delta_v = Kp*derivative(iterNum);
        v_dbi(iterNum+1) = v_dbi(iterNum) + delta_v;
    end
end
v_dbi = v_dbi(1:end-1); %clip off the last value because it's not real

figure(1);
plot(1:iterNum, maxClOverCdIterate, '-bo');
xlabel("Iteration Number"); ylabel("CL/CD"); 
title("Derivative-Based Iteration: Cl/CD vs. Iteration Number");
figure(2);
plot(1:iterNum, v_dbi, '-bo');
xlabel("Iteration Number"); ylabel("Velocity (m/s)"); 
title("Derivative-Based Iteration: Velocity vs. Iteration Number");

%update load factor
nIter = sqrt((v_dbi(iterNum)^2/radius/g)^2 + 1);
%update stall speed
stallSpeedIter = v_dbi(iterNum)/sqrt(maxLoadFactorStall);

maxClOverCd_sol = maxClOverCdIterate(iterNum);
vel_sol = v_dbi(iterNum);
best_cl_sol = bestcl(iterNum);
best_S_sol = bestS(iterNum);

%Update weight and lift
weight = getWeight(best_S_sol);
lift = weight;

wing_loading = lift/bestS(iterNum);
best_chord_sol = bestchord(iterNum);
best_AR_sol = bestAR(iterNum);

[tail_area_h, tail_area_v, tail_boom_length] = find_tail_size(wingSpan, ...
    best_S_sol, best_chord_sol, C_HT, C_VT, density, viscosity, vel_sol);

disp("----------------------------");
disp("DBI SOLUTION");
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
% 
%% Iterative solution
function [maxClOverCd_sol, vel_sol, best_cl_sol, best_S_sol, wing_loading, best_chord_sol, best_AR_sol, nIter, stallSpeedIter] ...
    = iter_sol(ARMaster, S, SMaster, cTail, density, g, lenFuse, lenNose, lenTailBoom, lift, maxClCruise, maxLoadFactorStall, maxLoadFactorTurns, radius, sFuse, sNose, sTail, sTailBoom, sweepAngle, taperRatio, velocity, viscosity, wingSpan, C_HT, C_VT) 
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

        %get the zero-lift drag coeff for this S and velocity
        cd0(i) = getZeroLiftDrag(density, viscosity, velocity(j), ...
                    S(i), chord(i), sFuse,lenFuse, sNose,lenNose,...
                    sTail,cTail, sTailBoom,lenTailBoom);

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
    best_S_sol, best_chord_sol, C_HT, C_VT, density, viscosity, vel_sol);

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

%% Analytic Solution
%Note: this is old as of November 2020. Run at your own risk.
function [maxClOverCd_sol, vel_sol, best_cl_sol, best_S_sol, wing_loading, best_chord_sol, best_AR_sol, nIter, stallSpeedIter] ...
    = anal_sol(cTail, density, g, lenFuse, lenNose, lenTailBoom, lift, maxLoadFactorStall, maxLoadFactorTurns, radius, sFuse, sNose, sTail, sTailBoom, sweepAngle, taperRatio, velocity, viscosity, wingSpan, C_HT, C_VT) 

%Iterate until S guess is close to S calc, and eGuess is close to e calc
%When comparing to the iterative solution, it appears that my assumption
%that finding the configuration with an optimal velocity of our desired
%velocity does not true for any velocity except the velocity which is the
%ideal velocity for the configuration. This can be seen in figure 6.

%%Preallocation for velocity (only iterates over velocity)
clOverCd = velocity;
cl = velocity;
cd0 = velocity;
cdi = velocity;
cd = velocity;
dragAna = velocity;
zeroLiftDragAna = velocity;
inducedDragAna = velocity;
s = velocity;
for i = 1:length(velocity)
    %%Declaring that we have these variables
    hasLift = 1;
    hasVel = 1;
    hasCd0 = 1;

    %want to find ideal S using minDragEq (source in function)
    sWingGuess = 0.8;    %initial guess based on Fiber One (old platform)
    s(i) = sWingGuess;
    sWingGuess = 0;     %so that is runs the loop
    infLoopCount = 0;   %counter to prevent infinite loop
    lastIter = false;   %signifies last iteration if inf loop triggered
    while(abs(sWingGuess - s(i)) > 10^-8)
        %Preparing for the getZeroLiftDrag function by stating we have the wing
        %area
        infLoopCount = infLoopCount + 1;
        if infLoopCount > 1000
            s(i) = (sWingGuess + s(i))/2;
            lastIter = true;
        else
            sWingGuess = s(i);
        end
        wingChord = s(i)/wingSpan;
        %get zero-lift drag coeff for this S
        cd0(i) = getZeroLiftDrag(density, viscosity, velocity(i), ...
        sWingGuess,wingChord, sFuse,lenFuse, sNose,lenNose, sTail,cTail, ...
        sTailBoom,lenTailBoom);

        %We want to solve for the wing area, so we state that we don't know it
        s(i) = -1;     %-1 means we don't have the information

        %Solve for ideal S and Cl for this velocity, wing span, cd0, and weight 
        [s(i), velocity(i), cl(i), k(i)] = minDragEq(velocity(i), density,...
            s(i), cd0(i), wingSpan, lift, taperRatio, sweepAngle);
        %this will return an s, but the cd0 was calculated based on a guess for
        %S, so this s becomes the new guess, and we iterate until the guess
        %matches the return. (causes oscillations sometimes on unlucky inputs)
        if lastIter
            break;
        end
    end

    %%Get characteristics
    %%Get cdi so we can later get cd
    cdi(i) = k(i)*cl(i)^2;    %cdi should equal cd0
    %%Get cd so we can later get cl/cd
    cd(i) = cdi(i) + cd0(i);
    %%Get drags for plotting
    dragAna(i) = 1/2*density*velocity(i)^2*s(i)*cd(i);
    zeroLiftDragAna(i) = 1/2*density*velocity(i)^2*s(i)*cd0(i);
    inducedDragAna(i) = 1/2*density*velocity(i)^2*s(i)*cdi(i);

    %%Get cl/cd to report performance
    clOverCd(i) = cl(i)/cd(i);
end

%%Plots for visualization
%%Analytic Results Plots
figure(201)
plot(velocity, dragAna, '-m', ...
    velocity, zeroLiftDragAna, '-k', velocity, inducedDragAna, '--c');
xlabel("Velocity (m/s)"); ylabel("Drag (N)"); 
title("Analytic: Minimum Drag vs. Velocity");
legend("total drag","zero lift drag","induced drag");
figure(202)
plot(velocity, clOverCd, '-b')
xlabel("Velocity (m/s)"); ylabel("CL/CD"); 
title("Analytic: Maximum CL/CD vs. Velocity");
%%Comparison plots between Iterative and Analytic
% figure(203)
% plot(velocity, bestDrag, '-g', ...
%     velocity, zeroLiftDrag, '-b', velocity, inducedDrag, '-r',...
%     velocity, dragAna, '--m', ...
%     velocity, zeroLiftDragAna, '-k', velocity, inducedDragAna, '--c');
% xlabel("Velocity (m/s)"); ylabel("Drag (N)"); 
% title("Iterative vs. Analytic: Minimum Drag vs. Velocity");
% legend("Iterative: total drag", "Iterative: zero lift drag", ...
%     "Iterative: induced drag",...
%     "Analytic: total drag", "Analytic: zero lift drag", ...
%     "Analytic: induced drag");
% figure(204)
% plot(velocity, maxClOverCdIterate, '-b', velocity, clOverCd, '--g')
% xlabel("Velocity (m/s)"); ylabel("CL/CD"); 
% title("Iterative vs. Analytic: Maximum CL/CD vs. Velocity");
% legend("Iterative", "Analytic");

%%Setting quantites to the most efficient
index = find(clOverCd == max(clOverCd));  %index of max efficiency
clOverCd = max(clOverCd);
velocityAna = velocity(index);
cl = cl(index);
cd = cd(index);
cd0 = cd0(index);
cdi = cdi(index);
s = s(index);
k = k(index);

chord = s/wingSpan;     %%Get chord to check feasibility
AR = wingSpan^2/s;  %%Get aspect ratio to check feasibility
nAna = sqrt((velocityAna^2/radius/g)^2 + 1); %%Update load factor for turns
stallSpeedAna = velocityAna/sqrt(maxLoadFactorStall);  %%Update stall speed

%Report characteristics
maxClOverCd_sol = clOverCd;
vel_sol = velocityAna;
best_cl_sol = cl;
best_S_sol = s;
wing_loading = lift/s;
best_chord_sol = chord;
best_AR_sol = AR;
nIter = nAna;
stallSpeedIter = stallSpeedAna;
disp("----------------------------");
disp("ANALYTIC SOLUTION");
disp("----------------------------");
disp("Max cl/cd: " + maxClOverCd_sol);
disp("Velocity: " + vel_sol);
disp("cl: " + best_cl_sol);
disp("S: " + best_S_sol);
disp("Wing loading: " + wing_loading);
disp("Chord: " + best_chord_sol);
disp("Aspect ratio: " + best_AR_sol);
disp("Load factor on turns: " + nIter);
if nIter > maxLoadFactorTurns + 0.001
    disp("Invalid solution! Load factor constraint not satisfied.");
end
disp("Stall speed: " + stallSpeedIter);
% %%Checks
% %%This method sets cd0 = cdi
% if abs(cdi - cd0) > 1e-5
%     disp("cd0 is not equal to cdi. Something went wrong");
% end
% %%This method should return a cl/cd of 1/(2*sqrt(cd0*k))
% expectedClOverCd = (1/(2*sqrt(cd0*k)));
% if abs(clOverCd - expectedClOverCd) > 1e-5
%     disp("cl/cd is not equal to expected value. Something went wrong.");
% end
% %%Lift should equal weight for level flight
% lift = 1/2*density*velocityAna^2*s*cl;
% if abs(lift - weight) > 1e-3
%     disp("lift is not equal to weight. Something went wrong.");
% end
end