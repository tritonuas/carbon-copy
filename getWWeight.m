function[weight] = getWWeight(area,length)
% Wing Area vs Weight
    %Assumptions Made:
        %Wing modeled as two flat plates; curvature not accounted for
        %Weight of epoxy and carbon equal
        %Two layers of carbon, top and bottom
        %Foam spar crossection modeled as .25" by 1"
% Assummed Weights:    
    %Divinycell
        d = 80; %oz/ft^3
              % https://www.fibreglast.com/product/vinyl-foam-5-lb-density/Foam  
    %Layup
        s = 2/9; %oz/ft^2
                % http://www.cstsales.com/carbon_prepreg.html 
% Calculation 
    %areastr = 'Wing Area for one wing, square feet ';
    %area = input(areastr); 
    
    %lengthstr = 'Wing length, tip to fuse, feet ';
    %length = input(lengthstr);
    
    ws = area*s*4;
    wb = .25*1/12*length*d;
    
    W = ws+wb;
    %disp('the weight in ounces is')
    weight = W
    
    
    
   
    
    