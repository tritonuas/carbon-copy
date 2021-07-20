density = 1.225;
viscosity = 1.789*10^-5;
velocity = 20;
sWing = 0.81547166;
saWing = 2*0.9576*sWing;
cWing = 0.22098;
saFuse = 0.9082*0.64;
saNose = 0;
lenFuse = 1.2;
saTail = 2*(0.082369 + 0.14259002); %% m^2
cVS = 0.300799;
cHS = 0.161552;
lenNose = -1;
saTailBoom = 2*pi*(0.0762/2);
lenTailBoom = 1;
saHS = 0.18932187;
saVS = 0.29549039;

%% Fudge Factors

wingFudgeFactor = 2.222293569;
fuseFudgeFactor = 1.954744348;
tailboomFudgeFactor = 0.4573414355;
tailFudgeFactor = 3.19497336;

%% cd0 calculations as a sum of the cd0s of all components. 

cd0 = getZeroLiftDrag(density, viscosity, velocity, ...
sWing, saWing,cWing, saFuse, lenFuse, saNose,lenNose,  ...
saHS,cHS, saVS,cVS,saTailBoom, lenTailBoom);

saWingcd0 = getZeroLiftDrag(density, viscosity, velocity, ...
sWing, saWing,cWing, 0, lenFuse, 0,lenNose,  ...
saHS,cHS, saVS,cVS,0, lenTailBoom);

saFusecd0 = getZeroLiftDrag(density, viscosity, velocity, ...
sWing, 0,cWing, saFuse, lenFuse, 0,lenNose,  ...
saHS,cHS, saVS,cVS,0, lenTailBoom);

saNosecd0 = getZeroLiftDrag(density, viscosity, velocity, ...
sWing, 0,cWing, 0, lenFuse, saNose,lenNose,  ...
saHS,cHS, saVS,cVS,0, lenTailBoom);

saTailcd0 = getZeroLiftDrag(density, viscosity, velocity, ...
sWing, 0,cWing, 0, lenFuse, 0,lenNose,  ...
saHS,cHS, saVS,cVS,0,lenTailBoom);

saTailboomcd0 = getZeroLiftDrag(density, viscosity, velocity, ...
sWing, 0,cWing, 0, lenFuse, 0,lenNose, ...
saHS,cHS, saVS,cVS,saTailBoom, lenTailBoom);

saWingcd0 = wingFudgeFactor*saWingcd0
saFusecd0 = fuseFudgeFactor*saFusecd0
saTailboomcd0 = tailboomFudgeFactor*saTailboomcd0
saTailcd0 = tailFudgeFactor*saTailcd0


saWingcd0+saFusecd0+saTailboomcd0+saTailcd0