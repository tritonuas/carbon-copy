function [Ix,Iy] = MomentofI(b1, b2, b3, h1, h2, h3)
%Triton UAS
%Moment of Inertia for an 'I' beam comprised of 3 rectangular sections
%Rectangles will be numbered from bottom to top
%Global y taken from base of material
%Global X is taken along the centerline of material
%User needs to specify: b1,b2,b3 and h1,h2,h3
%We are assuming rectangular sections so any additional section via fillet
%should be accounted for additional from this calculation
%Rectangle 1
A1=b1*h1;%area of rectangular section
Ix1=(1/12)*b1*h1^3;%follows as individual section moment of interia, X
Iy1=(1/12)*h1*b1^3;%follows as individual section moment of interia, Y
C1=h1/2; %follows as centroid of section 
D1=C1+(h2/2);
PA1=A1*D1^2; %Perpendicular Axis Therom
%Rectangle 2
A2=b2*h2; %area of rectangular section
Ix2=(1/12)*b2*h2^3;
Iy2=(1/12)*h2*b2^3;
%Rectangle 3
A3=b3*h3; %area of rectangular section
Ix3=(1/12)*b3*h3^3;
Iy3=(1/12)*h3*b3^3;
C3=(h3/2);
D3=C3+(h2/2);
PA3=A3*D3^2; %Perpendicular Axis Therom
%Global Ix, Iy
Ix=Ix2+Ix1+PA1+Ix3+PA3 %Calculation for Global Ix
Iy=Iy1+Iy2+Iy3 %Calculation for Global Iy
end