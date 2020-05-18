%%This is the main script (executable) for the plane design program
%%As of now, this takes in inputs of velocity, wingspan, and weight, and
%%calculates the ideal lift coefficient and wing reference area to maximize
%%cl/cd and therefore maximize range and minimize drag for cruise
%%conditions.
%%There are currently 2 approaches to this. One is an analyic method as
%%described in Anderson's 155A lecture slides of setting cdi=cd0, and this
%%is done using the minDragEq method. The other method is a purely 
%%iterative method for the purpose of checking. Currently, the iterative
%%process also can iterate over velocity to try to find an ideal velocity
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
format compact
format shortg;
addpath("Utils")

%% Constants
g = 9.81;                   %sea level Earth gravity
density = 1.225;            %sea level air density
viscosity = 1.789*10^-5;    %sea level air viscosity

%% Configuration
%subsonic plane
sweepAngle = 0;     %In the future, I can make this a function of velocity
taperRatio = getTaper(sweepAngle);  %see function for source of eq

%% Inputs
% I might make this whole program a function in the future
weight = 135;
lift = weight;      %want to design this for cruise conditions
wingSpan = 3.65;
velocity = 20;

%Fuselage Requirements
sFuse = 1;
lenFuse = 1.2;

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
%none of this is currently being used in the program, it's more for
%reference
radius = 20;    %want this radius
n = sqrt((velocity^2/radius/g)^2 + 1);%load factor needed for desired turns
if n > 2    %max load factor of 2
    n = 2;
    velocity = sqrt(radius*g*sqrt(n^2-1));  %desired velocity for turns
end

%% Iterative solution
%declaring that I have these variables for use of the lift equation
hasLift = 1;
hasDensity = 1;
hasVel = 1;

% S = linspace(.1,1,1000); %iterate over these S
S = 0.23;
chord = S/wingSpan;     %for fixed wing span
hasS = 1;       %declaring that I have S and want to calculate for Cl
hasCl = 0;      %This is in preperation for the lift equation
AR = getAR(wingSpan, S);  %Getting ARs for induced drag calculations


%%preallocation for inner loop (S)
cl = S;
cd = S;
cd0 = S;
cdi = S;
clOverCd = S;
drag = S;
liftGen = S;

% velocity = 10:.01:40; %iterate over these velocity
velocity = 20.675;
%%preallocation for outer loop (velocity)
maxClOverCdIterate = velocity;
bestS = velocity;
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
    
    for i = 1:length(S) %for every wing area...
        
        %get the lift coefficient necessary for this S and velocity
        [lift, density, velocity(j), S(i), cl(i)] = liftEq(hasLift,lift,...
           hasDensity,density, hasVel,velocity(j), hasS,S(i), hasCl,cl(1));

        %get the zero-lift drag coeff for this S and velocity
        cd0(i) = getZeroLiftDrag(density, viscosity, velocity(j), ...
                    S(i), chord(i), sFuse,lenFuse, sNose,lenNose,...
                    sTail,cTail, sTailBoom,lenTailBoom);

        %get induced drag so we can later get cd
        k = getK(AR(i), taperRatio, sweepAngle);
        cdi(i) = cl(i)^2*k;
        %get cd so we can later get cl/cd
        cd(i) = getDragCoeff(cd0(i),cdi(i));
        %get cl/cd for performance comparison
        clOverCd(i) = cl(i)/cd(i);
        %get max cl/cd for best S for each velocity
        if clOverCd(i) >= maxClOverCd
            maxClOverCd = clOverCd(i);
            maxIndex = i;
        end
        drag(i) = dragEq(0,0, hasDensity,density, hasVel,velocity(j),...
            hasS,S(i), 1,cd(i));
        
    end
    
    %get characteristics of each local best performance
    maxClOverCdIterate(j) = maxClOverCd;
    bestS(j) = S(maxIndex);
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

%%Plots for visualization
figure(1)
plot(velocity, bestDrag, '-k', ...
    velocity, zeroLiftDrag, '-r', velocity, inducedDrag, '-m')
xlabel("Velocity (m/s)"); ylabel("Drag (N)"); 
title("Minimum Drag vs. Velocity");
legend("total drag","zero lift drag","induced drag");
figure(2)
plot(velocity, maxClOverCdIterate, '-b')
xlabel("Velocity (m/s)"); ylabel("CL/CD"); 
title("Maximum CL/CD vs. Velocity");

disp("----------------------------");
disp("ITERATIVE SOLUTION");
disp("----------------------------");
disp("Max cl/cd: " + maxOfMaxClOverCd);
disp("Best Velocity: " + velocity(maxOfMaxIndex));
disp("Best cl: " + bestcl(maxOfMaxIndex));
disp("Best S: " + bestS(maxOfMaxIndex));

%% Analytic Solution
%Maybe I should iterate until the sGuess = the sCalculated
%maybeI could also do that with the e inside the function
%Rn for our purposes, its giving me an s of about 0.3588m^2, which gives a
%wingChord of about 0.0983m (~0.322 feet) and an aspect ratio of about 37.1
%STILL NEEDS WORK. Wildly sensitive solutions. For whatever reason, it spits
%out a reasonable solution for v = 18.

%%INPUTS:
velocity = 20.765;

%%Declaring that we have these variables
hasLift = 1;
hasVel = 1;
hasCd0 = 1;

%want to find ideal S using minDragEq (source in there)
sWingGuess = 0.8;    %initial guess based on Fiber One (old platform)
s = sWingGuess;
sWingGuess = 0;     %so that is runs the loop
while(abs(sWingGuess - s) > 10^-7)
    
    %Preparing for the getZeroLiftDrag function by stating we have the wing
    %area, but do not have any of these other variables (guesses will be
    %given for variables we do not have from within the function)
    %-1 means we do not have the information
    sWingGuess = s;
    wingChord = s/wingSpan;
    sFuse = -1;
    lenFuse = -1;
    sNose = -1;
    lenNose = -1;
    sTail = -1;
    cTail = -1;
    sTailBoom = -1;
    lenTailBoom = -1;
    
    %get zero-lift drag coeff for this S
    cd0 = getZeroLiftDrag(density, viscosity, velocity, ...
    sWingGuess,wingChord, sFuse,lenFuse, sNose,lenNose, sTail,cTail, ...
    sTailBoom,lenTailBoom);

    %We want to solve for the wing area, so we state that we don't know it
    s = -1;     %-1 means we don't have the information
    
    %Solve for ideal S and Cl for this velocity, wing span, cd0, and weight 
    [s, velocity, cl, k] = minDragEq(velocity, density, s, cd0, ...
        wingSpan, lift, taperRatio, sweepAngle);
    %this will return an s, but the cd0 was calculated based on a guess for
    %S, so this s becomes the new guess, and we iterate until the guess
    %matches the return. (causes oscillations sometimes on unlucky inputs)
end

%%Get characteristics

%%Get cdi so we can later get cd
cdi = getInducedDragCoeff(k,cl);    %cdi should equal cd0
%%Get cd so we can later get cl/cd
cd = getDragCoeff(cd0, cdi);
%%Get cl/cd to report performance
clOverCd = cl/cd;

%%Get chord to check feasibility
chord = chordFromS(s, wingSpan);
%%Get aspect ratio to check feasibility
AR = wingSpan^2/s;

%%Report characteristics
disp("----------------------------");
disp("ANALYTIC SOLUTION");
disp("----------------------------");
disp("cl/cd: " + clOverCd);
disp("Velocity (input echo): " + velocity + " m/s");
disp("cl: " + cl);
disp("S: " + s + " m^2");
disp("Chord: " + chord + " m");
disp("Aspect Ratio: " + AR);

%%Checks
%%This method sets cd0 = cdi
if abs(cdi - cd0) > 1e-5
    disp("cd0 is not equal to cdi. Something went wrong");
end
%%This method should return a cl/cd of 1/(2*sqrt(cd0*k))
expectedClOverCd = (1/(2*sqrt(cd0*k)));
if abs(clOverCd - expectedClOverCd) > 1e-5
    disp("cl/cd is not equal to expected value. Something went wrong.");
end
%%Lift should equal weight for level flight
if abs(lift - weight) > 1e-5 
    disp("lift is not equal to weight. Something went wrong.");
end






