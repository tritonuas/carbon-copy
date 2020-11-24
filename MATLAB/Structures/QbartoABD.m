%%Finds the ABD matrices of a laminate
%%Author: Andrew Fletcher
%%Email: afletche@ucsd.edu
%@param Qbars the 3D matrix of Qbar matrices for each lamina
%@param thicknesses the vector of thicknesses for each lamina
%@return ABD the complete ABD matrix for the laminate
%@return A the A (extensional) matrix for the laminate
%@return B the B (in-plnae/out-of-plane coupling) matrix for the laminate
%@return D the D (bending stiffness) matrix for the laminate
%@return Ex the effective moduli, Ex, if symmetric layup with in plane load
%@return Ey " "                 , Ey, " "
%@return G12 " "                , Gxy, " "
%@return Nu12 " "               , Nuxy, " "
function [ABD, A, B, D, Ex, Ey, Gxy, Nuxy] = QbartoABD(Qbars, thicknesses)

%find height of the middle of each lamina
total_height = sum(thicknesses);
mid_plane_height = total_height/2;
z = cumsum(thicknesses) - mid_plane_height - thicknesses./2;

%Preallocation
A = zeros(length(Qbars(:,1,1)), length(Qbars(1,:,1)));
B = zeros(length(Qbars(:,1,1)), length(Qbars(1,:,1)));
D = zeros(length(Qbars(:,1,1)), length(Qbars(1,:,1)));

for k = 1:length(Qbars(1,1,:))
    %finding lower and upper z's
    z_top = z(k) + thicknesses(k)/2;
    z_bot = z(k) - thicknesses(k)/2;
    
    %Integral definition of A
    %A = A + Qbars(:,:,k)*(z_top - z_bot)
    A = A + Qbars(:,:,k)*thicknesses(k);     %slightly more efficient
    %Integral definition of B
    %B = B + Qbars(:,:,k)*1/2*(z_top^2 - z_bot^2);
    B = B + Qbars(:,:,k)*(z_top^2 - z_bot^2)/2; %slightly more efficient
    %Integral definition of D
    %D = D + Qbars(:,:,k)*1/3*(z_top^3 - z_bot^3);
    D = D + Qbars(:,:,k)*(z_top^3 - z_bot^3)/3;  %slightly more efficient
end
    
ABD = [A B; B D];
Ex = 1/(A(1,1)*total_height);
Ey = 1/(A(2,2)*total_height);
Gxy = 1/(A(3,3)*total_height);
Nuxy = -A(1,2)/A(1,1);