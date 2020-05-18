%%Author: Andrew Fletcher
%%Email: afletche@ucsd.edu
%%This program for finding the induced drag factor, K.
%%This is the K found in Cdi = kCl^2

%%This may seem ridiculous, but I am someone who makes many stupid errors,
%%and this should be a good way of reducing them. Especially since I seem
%%to do this computation fairly often.

%@param ar the aspect ratio, wingspan^2/wingArea
%@param taperRatio the taper ratio, chordTip/chordRoot
%@return k the induced drag factor
function [k] = getK(ar, taperRatio, sweepAngle)

e = getOswaldEff(taperRatio, ar, sweepAngle);
k = 1/(pi*e*ar);

end

%for finding the oswalds efficiency factor, e
%should maybe be taken with a grain of salt.
%If possible, I would recommend more research into this topic
%@param taperRatio the taperRatio (ct/cr)
%@param ar the aspect ratio (b^2/s)
%@param velocity the velocity of the aircraft
%@return e the oswalds efficiecy factor
function e = findOswaldEff(sweepAngle, ar)
%Taken from a mix of eq [1] and eq [21] from:
%https://www.fzt.haw-hamburg.de/pers/Scholz/OPerA/OPerA_PUB_DLRK_12-09-10.pdf
%seems to give pretty unrealistic values imo
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
%location being the main consequences of taper ratio).
%e = 1.78*(1 - 0.045*ar^(0.68)) - 0.64;

%According to the twice already cited paper, the optimal taper ratio is
%given by the equation
%sweepAngle = 0;      %main purpose of sweeping is to avoid drag from transonic
                    %and supersonic speeds. We are well below transonic
lambdaOpt = getTaper(sweepAngle); %sweep of quarter chord
%This gives an optimal taper ratio = 0.45 for no sweep (seems reasonable
%and the logic and result compare very closely to what I have read before)
%clarification: lambda = taper ratio

%the paper then goes on to give an equation for theoretical oswald factor
func = @(x) 0.0524*x^4 - 0.15*x^3 + 0.1659*x^2 - 0.0706*x + 0.0119;
%min of that function is x = 0.357
eTheoretical = 1/(1 + func(0.357)*ar);   %for min taper ratio
%actual equation is given in the form eTheo = 1/(1 + func(lam - delLam)*ar)
%where delLam is -0.357 + 0.45*exp(0.0375*sweep)

%then three correction factors are given, keF for fuselage influence, FeD0
%for zero-lift drag, and keM for mach number influence.
%approx values are given in TAB 3. of the paper:
%for general aviation:
keF = 0.971;    %table
keD0 = 0.804;   %table
keM = 1;        %for subsonic aircraft
e = eTheoretical*keF*keD0*keM;

%The paper goes on to find corrections for winglets, dihedral, and more
%unconventional wing configurations such as box wings and c-wings.
%I will leave them out of the calculations for now since we/I have not
%determined the dihedral yet, nor whether we will actually have winglets
%also because adding the dihedral would require me to tweak the method to
%ask for wingspan and dihedral angle when it would otherwise not be needed
%kWL is the penalization factor used to make the data fit real life. There
%is a lot of data scatter for specific examples, so I used the average
%which is pretty concerningly arbitrary.
%kWL = 2.83;
%keDihedral = (1 + 1/kWL*(1/cos(dihedralAngle - 1))^2
%e = e*Theoretical*KeF*KeD0*KeM*KeDihedral
end


