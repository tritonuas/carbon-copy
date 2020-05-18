%% General
format compact
addpath("Utils")

%% Iterative solution
%%This is for testing other files
g = 9.81;

hasLift = 1;
lift = 135;
wingSpan = 3.65;
hasDensity = 1;
density = 1.225;
viscosity = 1.789*10^-5;
radius = 30;
n = 1.55;
velocity = sqrt(radius*g*sqrt(n^2-1)) %calculating desired vel for turns
hasVel = 1;
S = linspace(.1,1,100);
hasS = 1;
hasCl = 0;
cl = S;
chord = S/wingSpan;
AR = (wingSpan^2)*ones(1,length(S))./S;

cd0 = S;
hasSWing = 1;
reWing = S;
hasSFuse = 0;
sFuse = 0;
reFuse = 0;
isLamFuse = 1;
hasSNose = 0;
sNose = 0;
reNose = 0;
isLamNose = 1;
hasSTail = 0;
sTail = 0;
reTail = 0;
isLamTail = 1;
hasSTailBoom = 0;
sTailBoom = 0;
reTailBoom = 0;
isLamTailBoom = 1;

cdi = S;
cd = S;

clOverCd = S;
maxClOverCd = 0;

drag = S;
liftGen = S;
velocity = 10:.1:50;
for j = 1:length(velocity)
    maxClOverCd = 0;
    cd0 = S;
    sFuse = -1;
    lenFuse = -1;
    sNose = -1;
    lenNose = -1;
    sTail = -1;
    cTail = -1;
    sTailBoom = -1;
    lenTailBoom = -1;
    S = linspace(.1,1,100);
    for i = 1:length(S)
        hasS = 1;
        hasCl = 0;
        [lift, density, velocity(j), s, cl(i)] = liftEq(hasLift,lift, ...
            hasDensity,density,hasVel,velocity(j), hasS,S(i), hasCl,cl(1));

        reWing(i) = getRe(density,velocity(j),chord(i),viscosity);
        if reWing(i) <= 250000      %update this for Rex over flat plate
            isLamWing = 1;
        else
            isLamWing = 0;
        end
        cd0(i) = getZeroLiftDrag(density, viscosity, velocity(j), ...
                    S(i),S(i)/wingSpan, sFuse,lenFuse, sNose,lenNose,...
                    sTail,cTail, sTailBoom,lenTailBoom);

        sweepAngle = 0;

        k = getK(AR(i),.45, sweepAngle);
        cdi(i) = cl(i)^2*k;

        cd(i) = cdi(i) + cd0(i);

        clOverCd(i) = cl(i)/cd(i);
        if clOverCd(i) >= maxClOverCd
            maxClOverCd = clOverCd(i);
            maxIndex = i;
        end

        drag(i) = cd(i)*S(i)*1/2*density*velocity(j)^2;
        liftGen(i) = cl(i)*S(i)*1/2*density*velocity(j)^2;
    end

    maxClOverCdIterate(j) = maxClOverCd;
    maxIndex;
    bestS(j) = S(maxIndex);
    bestcd0(j) = cd0(maxIndex);
    bestcdi(j) = cdi(maxIndex);
    bestcd(j) = cd(maxIndex);
    bestcl(j) = cl(maxIndex);
    bestDrag(j) = drag(maxIndex);
    bestLift(j) = liftGen(maxIndex);
end
figure(1)
plot(velocity, bestDrag)

figure(2)
plot(velocity, maxClOverCdIterate)

%% Testing minDragEq
%Maybe I should iterate until the sGuess = the sCalculated
%maybeI could also do that with the e inside the function
%Rn for our purposes, its giving me an s of about 0.3588m^2, which gives a
%wingChord of about 0.0983m (~0.322 feet) and an aspect ratio of about 37.1
%STILL NEEDS WORK. Wildly sensitive solutions. For whatever reason, it spits
%out a reasonable solution for v = 20, but I don't think that should be
%trusted.

%inputs
hasViscosity = 1;
viscosity = 1.789*10^-5;    %air at sea level in SI
hasDensity = 1;
density = 1.225;            %air at sea level in SI
hasLift = 1;
lift = 135;                 %our guess for the weight
hasVelocity = 1;
velocity = 20;              %our guess for a good flight velocity
%velocity will be bounded by turn radius, pictures, ability to finish
%competition on time, and feasibility in aerodynamics
hasWingSpan = 1;
wingSpan = 3.65;            %assuming we want the max wing span we can get
%wingspan is bounded by out ability to transport it (and theoretically
%runway size). We assume the max is the ideal since the added weight of
%structure is probably negligible compared to aerodynamic benefit thanks to
%composites

hasCd0 = 1;
sWingGuess = 1;    %will iterate
hasSWing = 1.1;
s = sWingGuess;
sWingGuess = 0;     %so that is runs the loop
while(abs(sWingGuess - s) > 10^-7)
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
    
    cd0 = getZeroLiftDrag(density, viscosity, velocity, ...
    sWingGuess,wingChord, sFuse,lenFuse, sNose,lenNose, sTail,cTail, ...
    sTailBoom,lenTailBoom);

    hasS = 0;
    s = 0;
    
    [s, velocity, k] = minDragEq(hasVelocity, velocity,...
        hasDensity, density, hasS, s, hasCd0, cd0, hasWingSpan, wingSpan, ...
        hasLift, lift);
end

s
chord = chordFromS(s, wingSpan)
AR = wingSpan^2/s

hasS = 1;
hasCl = 0;
cl = 0;
[lift, rho, velocity, s, cl] = liftEq(hasLift,lift, hasDensity,density,...
    hasVelocity,velocity, hasS,s, hasCl,cl);

cdi = getInducedDragCoeff(k,cl);    %got it such that cdi = cd0
cd = cdi + cd0
cdi
cd0
clOverCd = cl/cd
cl
drag = cd*1/2*density*velocity^2*s;
lift = cl*1/2*density*velocity^2*s;
