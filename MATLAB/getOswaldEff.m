%%Author: Andrew Fletcher
%%Email: afletche@ucsd.edu
%%This program for finding the Oswald Efficiency factor, e

%This program follows several of the methods listed out by a research
%paper that can be found here:
%https://www.fzt.haw-hamburg.de/pers/Scholz/OPerA/OPerA_PUB_DLRK_12-09-10.pdf
%I have reasonable trust for this source as they have seem to have done
%their research, and what logic I picked through seemed to make sense, but
%it is still only one source.
%Additionally, I was decently careful to only use equations that contained
%the same assumptions (or fewer assumptions) than I was making.
%I recommend doing more research into making sure the equation checks out,
%ESPECIALLY if you have a sweep angle ~= 0.

%@param taperRatio the taperRatio (ct/cr)
%@param ar the aspect ratio (b^2/s)
%@param sweepAngle the angle of the sweep of the wings at quarter chord
%@return e the oswalds efficiecy factor
function e = getOswaldEff(taperRatio, ar, sweepAngle)
%Taken from a mix of eq [1] and eq [21] from:
%https://www.fzt.haw-hamburg.de/pers/Scholz/OPerA/OPerA_PUB_DLRK_12-09-10.pdf
%seems to give pretty unrealistic values, rip
% sos = 344;  %speed of sound at sea level
% machN = velocity/sos
% beta = sqrt(1 - machN^2)
% sigma = (0.0015 + 0.014*(taperRatio - 0.4)^2)*(beta*ar-4.5);    %eq [21]
% e = 1/(1 + sigma)  %eq[1]

%From MAE 155A notes from Anderson (also found in):
%https://www.fzt.haw-hamburg.de/pers/Scholz/OPerA/OPerA_PUB_DLRK_12-09-10.pdf
%as eq 32 from the 21st author (Raymer). assumptions are for unswept wings
%and sweep < 30 degrees with ar > 2.27 (which gives e > 1)
%My only problem with this formula is it does not take into account
%taperRatio, which according to my prior reading, taperRatio seemed to be
%the a heavy factor (with this efficiency factor and the initial stall
%location being the main consequences of taper ratio). Additionally, the
%paper declares that this method, while empirically calculated, is only
%applicable in a limited design space, (whatever that means), and 'error
%quite large'.
%e = 1.78*(1 - 0.045*ar^(0.68)) - 0.64;


%the paper then goes on to give it's own equation for theoretical e
func = @(x) 0.0524*x^4 - 0.15*x^3 + 0.1659*x^2 - 0.0706*x + 0.0119;
%min of that function is at x = 0.357
%eTheoretical = 1/(1 + func(0.357)*ar)   %assuming sweep = 0
%actual equation is given in the form:  FISHY for any sweep ~= 0!!
taperRatioOpt = getTaper(sweepAngle);
deltaTaperRatio = -0.357 + taperRatioOpt;
eTheoretical = 1/(1 + func(taperRatio - deltaTaperRatio)*ar);

%then three correction factors are given, keF for fuselage influence, FeD0
%for zero-lift drag, and keM for mach number influence.
%approx values are given in TAB 3. of the paper:
keF = 0.971;    %table value for general aviation
keD0 = 0.804;   %table value for general aviation
keM = 1;        %for subsonic aircraft
% e = eTheoretical*keF*keD0*keM;
e = eTheoretical;

%The paper goes on to find corrections for winglets, dihedral, and more
%unconventional wing configurations such as box wings and c-wings.
%I will leave them out of the calculations for now since we/I have not
%determined the dihedral yet, nor whether we will actually have winglets.
%also, adding the dihedral would require me to tweak the method to
%ask for dihedral angle when it would otherwise not be needed
%kWL is the penalization factor used to make the data fit empirical 
%results. There is a lot of data scatter for specific examples, so I used
%the average which is pretty concerningly arbitrary.
%kWL = 2.83;
%keDihedral = (1 + 1/kWL*(1/cos(dihedralAngle - 1))^2
%e = e*Theoretical*KeF*KeD0*KeM*KeDihedral
end