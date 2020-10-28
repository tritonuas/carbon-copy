%%Gets the weight of the plane, based on an unknown weight of the plane
%%body + wing weight, which depends on area S.
%%Weight = area*(s*4+d*t*2)*g, where g = 9.81, area = wing area
%%s = layup, (d, t) = density / thickness of divinycell
function[weight] = getWeight(area,length)

g = 9.81; %acc due to gravity
% Wing Area vs Weight
    %Assumptions Made:
        %Wing modeled as two flat plates; curvature not accounted for
        %Weight of epoxy and carbon equal
        %Two layers of carbon, top and bottom
        %Foam spar crossection modeled as .25" by 1"
% Assummed Weights:    
    %fiberOne wing weight: around 2 pounds (0.907kg)
    % assumed weight is 135N (13.76kg)
    wPlane = (135/g) - 0.907;


    %Divinycell
        %d = 80; oz/ft^3
        d = 80.1; % (1 oz/ft^3 = 1.001kg/m^3)
        %thickness = .25"
        t = 0.00635; % in m
              % https://www.fibreglast.com/product/vinyl-foam-5-lb-density/Foam  
    %Layup
        %s = 2/9; oz/ft^2
        s = 0.0678; % (1 oz/ft^2 = 0.305 kg/m^2)
                % http://www.cstsales.com/carbon_prepreg.html 
% Calculation 
    %areastr = 'Wing Area for one wing, square feet ';
    %area = input(areastr); 
    
    %lengthstr = 'Wing length, tip to fuse, feet ';
    %length = input(lengthstr);
    
    ws = area*s*4; % weight of skin
   % wb = .25*1/12*length*d; % weight of rib, .25 and 1/12 in ft
    wb = area*d*t; % weight of full divinycell wing
   m = 2*(ws+wb) + wPlane;
   
   %weight = area*(s*4+d*t*2); ??? forgot what this was
    %disp('the weight in N is')
   
    weight = m*g;
    