%%This is for testing other files
% Units: meters, seconds, Kg, N

hasLift = 1;
lift = 266.984447; % 60 lb = 27.21554 kg, 27.21554 kg * 9.81 m/s^2 = 266.984447 N
wingSpan = 4.1148; % 2 ft * 2 + 114 inches = 13.5 ft = 4.1148 m
hasRho = 1;
rho = 1.225;
viscosity = 1.789*10^-5;
hasVel = 1;
velocity = 40.2336; % 132 ft/s
hasS = 1;
S = linspace(.1,3,100);
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
for i = 1:length(S)
    [lift, rho, velocity, s, cl(i)] = liftEq(hasLift,lift, hasRho,rho,...
    hasVel,velocity, hasS,S(i), hasCl,cl(1));
    
    reWing(i) = getRe(rho,velocity,chord(i),viscosity);
    if reWing(i) <= 250000
        isLamWing = 1;
    else
        isLamWing = 0;
    end
    cd0(i) = getZeroLiftDrag(hasSWing,S(i),reWing(i),...
    hasSFuse,sFuse,reFuse, hasSNose,sNose,reNose, hasSTail,sTail,reTail,...
    hasSTailBoom,sTailBoom,reTailBoom, isLamWing, isLamFuse,...
    isLamNose, isLamTail, isLamTailBoom);

    sweepAngle = 0;

    k = getK(AR(i),.45, sweepAngle);
    cdi(i) = cl(i)^2*k;
    
    cd(i) = cdi(i) + cd0(i);
    
    clOverCd(i) = cl(i)/cd(i);
    if clOverCd(i) >= maxClOverCd
        maxClOverCd = clOverCd(i)
        maxIndex = i
    end
    
    drag(i) = cd(i)*S(i)*1/2*rho*velocity^2;
    liftGen(i) = cl(i)*S(i)*1/2*rho*velocity^2;
end

%% Testing minDragEq
%Maybe I should iterate until the sGuess = the sCalculated
%maybeI could also do that with the e inside the function
%Rn for our purposes, its giving me an s of about 1.5315, which gives a
%wingChord of about 0.4196m (~1.5 feet) and an aspect ratio of about 8.7

viscosity = 1.789*10^-5;
density = 1.225;

hasVelocity = 1;
velocity = 40.2336;
hasDensity = 1;
density = 1.225;
hasS = 0;
s = 0;

hasCd0 = 0;
hasSWing = 1;
sWing = 0.2955;
wingSpan = 4.1148;
wingChord = sWing/wingSpan;
reWing = getRe(density, velocity, wingChord, viscosity);
hasSFuse = 0;
sFuse = 0;
reFuse = 0;
hasSNose = 0;
sNose = 0;
reNose = 0;
hasSTailBoom = 0;
sTailBoom = 0;
reTailBoom = 0;
isLamWing = 0;
isLamFuse = 0;
isLamNose = 0;
isLamTail = 0;
isLamTailBoom = 0;
cd0 = getZeroLiftDrag(hasSWing,sWing,reWing,...
    hasSFuse,sFuse,reFuse, hasSNose,sNose,reNose, hasSTail,sTail,reTail,...
    hasSTailBoom,sTailBoom,reTailBoom, isLamWing, isLamFuse,...
    isLamNose, isLamTail, isLamTailBoom);
hasWingSpan = 1;
wingSpan = 4.1148;
hasLift = 1;
lift = 266.984447

[s, velocity] = minDragEq(hasVelocity, velocity,...
    hasDensity, density, hasS, s, hasCd0, cd0, hasWingSpan, wingSpan, ...
    hasLift, lift)

C_L = 2*lift/(s*rho*velocity^2)
c = s/wingSpan;
disp(c * 39.37008); % Convert from m to inches
disp("in");


% hasLift = 1;
% hasRho = 1;
% rho = 1.225;
% hasVel = 0;
% velocity = 0;
% hasS = 1;
% hasCl = 1;
% cl = 1.5
% [lift, rho, velocity, s, cl] = liftEq(hasLift,lift, hasRho,rho,...
%     hasVel,velocity, hasS,s, hasCl,cl)
