%%Author: Andrew Fletcher
%%Email: afletche@ucsd.edu
%%This program for finding the zero-lift component of the drag, often
%%denoted as cd0.

%SOURCE for equations: Professor Anderson (MAE 155A lecture slides),
%Munson, Young, Okishi Fluid Mechanics 8th edition, and MAE 101C Professor 
%Bahadori (coming from Basic Heat Transfer, Mills 3rd edition) all use the 
%same equations. However, Prfoessor Hwang, citing Raymer Aviation Design
%uses a different process.
%ASSUMPTIONS: The laminar solution is using the Blassius solution, which 
%assumes a smooth surface and laminar flow. The turbulent solution is 
%also empirical and assumes a turbulent flow and completely smooth surface. 
%OTHER NOTES on equation(s): they are empirically derived equations
%CONCERNS: This equation only seems to account for the skin friction drag,
%and not any other factors such as cross-sectional area drag. Im not too
%well read on this, but for my current understanding, that is a slight
%concern. Also, equation does not account for fuse tapering toward tail,
%nor does it account for tail boom? Tail boom portion is my addition. I
%think it is right, but Im not 100% sure.
%CONCERNS RESOLVED: Usually, only the skin friction is taken into account 
%for a thin or streamlined body since the cross sectional area is so low.
%Anyways, the cross sectional area drag is actually just the pressure drag,
%which is taken into account by induced drag (D0 is for skin frction, Di is
%for pressure drag). As for the equations, see Munson, Young, Okishi Fluids
%Mechanics 8th edition, section 9.3. These are empirical formulas, but I 
%believe there are no known alternatives, and I believe these are fairly
%well established.
%More research recommended?: low priority, but yes

%@param density the density of the air
%@param viscosity the viscosity of the air
%@param velocity the velocity to the flow relative to body
%@param sWing the wetted area of the wings
%@param cWing the average chord of the wings
%@param sFuse the wetted area of the fuselage
%@param lenFuse the length of the fuselage
%@param sNose the wetted area of the nose
%@param lenNose the length of the nose
%@param sTail the wetted area of the tail
%@param cTail the average chord of the tail
%@param sTailBoom the wetted area of the tail boom
%@param lenTailBoom the length of the tail boom
%@return cD0 the parasite drag coefficient (zero-lift drag coefficient)
function [cd0] = getZeroLiftDrag(density, viscosity, velocity, ...
    sWing,cWing, sFuse,lenFuse, sNose,lenNose, sTail,cTail, ...
    sTailBoom,lenTailBoom)
%guesses for a size plane Triton UAS would be using
%only used if values are not given in parameters (input as -1)
%Guesses are given based off of my roughish guesses for Fiber One
if sWing == -1
    sWing = 2;         %all in m^2
end
if cWing == -1
    cWing = 0.6;
end

if sFuse == -1
    sFuse = 0.4625;
end
if lenFuse == -1
    lenFuse = 0.8;   %m
end

if sNose == -1
    sNose = 0.0325;
end
if lenNose == -1
    lenNose = 0.1;   %m
end

if sTail == -1
    sTail = 0.42;
end
if cTail == -1
    cTail = 0.2;   %m
end

if sTailBoom == -1
    sTailBoom = 0.0997;
end
if lenTailBoom == -1
    lenTailBoom = 0.9;   %m
end

if velocity == -1
    velocity = 20;     %m/s
end
if density == -1
    density = 1.225;   %kg/(m^3)
end
if viscosity == -1
    viscosity = 1.789*10^-5;   %kg/(m*s);
end

cfWing = getCf(velocity, density, viscosity, cWing);
cfFuse = getCf(velocity, density, viscosity, lenFuse);
cfNose = getCf(velocity, density, viscosity, lenNose);
cfTail = getCf(velocity, density, viscosity, cTail);
cfTailBoom = getCf(velocity, density, viscosity, lenTailBoom);

areaRatioWing = sWing/sWing;
areaRatioFuse = sFuse/sWing;
areaRatioNose = sNose/sWing;
areaRatioTail = sTail/sWing;
areaRatioTailBoom = sTailBoom/sWing;

cd0 = 1.25*(cfWing*areaRatioWing + cfFuse*areaRatioFuse + ...
    cfNose*2/sqrt(3)*areaRatioNose + cfTail*areaRatioTail + ...
    cfTailBoom*areaRatioTailBoom);

%@param re Reynold's number
%@param isLaminar   boolean for whether the surface has laminar flow or
%turbulent
function cf = getCf(velocity, density, viscosity, chord)
Re = density*velocity*chord/viscosity;

RE_CRIT_TRANS = 1*10^5; %currently unused because no transition flow
% Original
% RE_CRIT_TURB = 2*10^5;
% Hwang's Method:
k = 0.7*10^-5;
RE_CRIT_TURB = 38.21*(chord/k)^1.053;

xCritTrans = RE_CRIT_TRANS*viscosity/(velocity*density);
xCritTurb = RE_CRIT_TURB*viscosity/(velocity*density);

lamFunc = @(x)(1.328/sqrt(density*velocity*x/viscosity));
transFunc = @(x)(0.455/(log10(density*velocity*x/viscosity)^2.58) ...
    - 1700/(density*velocity*x/viscosity));
turbFunc = @(x)(0.455/(log10(density*velocity*x/viscosity)^2.58));

turbFuncMills = @(x)(lamFunc(xCritTurb)*RE_CRIT_TURB/x + ...
    0.523/log(0.006*x)^2*(1-RE_CRIT_TURB/x));
turbFuncMillsLessThanE7 = @(x)(lamFunc(xCritTurb)*RE_CRIT_TURB/x + ...
    0.0174*x^(-1/2)*(1-(RE_CRIT_TURB/x)^(4/5)));

% %Mills Method
% if Re < RE_CRIT_TURB
%     cf = lamFunc(chord);
% % elseif Re < 10^7
% %     cf = turbFuncMillsLessThamE7(Re);
% else
%     cf = turbFuncMills(Re);
% end

fullLamComponent = lamFunc(xCritTrans);
% fullTransComponent = transFunc(xCritTurb) - transFunc(xCritTrans);
if chord > xCritTurb
    cf = (fullLamComponent*xCritTurb + ...
          turbFunc(chord)*(chord-xCritTurb))...
        /chord;
% elseif chord > xCritTrans
%     cf = (fullLamComponent*xCritTrans + ...
%         transFunc(chord)*(chord-xCritTrans))...
%         /chord;
else
    cf = lamFunc(chord);
end
