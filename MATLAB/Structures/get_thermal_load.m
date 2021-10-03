%%Takes in inputs of material and laminate properties and calculates the
%%thermal loading.
%%Author: Andrew Fletcher
%%Email: afletche@ucsd.edu
%
%%@param thetas the angle of each ply in the laminate from bottom to top
%%@param rad_or_deg the string 'rad' or 'deg' specifying units of thetas
%%@param mat_props a vector of material properties defined as [E11 E22 G12 Nu12]
%%@param thicknesses a vector of the thickness of each ply
%%@param delta_T the change in temperature that causes the thermal loading
%%@param cte_mat_local the matrix of coefficients of thermal expansion
%%@return thermal_loading the applied load from temperature changes
function thermal_loading = get_thermal_load(thetas, rad_or_deg,...
    mat_props, thicknesses, delta_T, cte_mat_local)

%find height of the middle of each lamina
total_height = sum(thicknesses);
mid_plane_height = total_height/2;
z = cumsum(thicknesses) - mid_plane_height - thicknesses./2;

%find Q_bars
for k = 1:length(thetas)
    T_sigma = get_T_planar(thetas(k), 'stress', 'deg');
    cte_mat_global(:,k) = T_sigma'...
        *cte_mat_local(:,k);
    S(:,:,k) = get_S_planar_ortho(mat_props(:,k));
    Q(:,:,k) = inv(S(:,:,k));
    Q_bars(:,:,k) = get_Q_bar_planar(Q(:,:,k), thetas(k), rad_or_deg);
end

%initialize
N_thermal = 0;
M_thermal = 0;

for k = 1:length(Q_bars(1,1,:))
    %finding lower and upper z's
    z_top = z(k) + thicknesses(k)/2;
    z_bot = z(k) - thicknesses(k)/2;
    
    N_thermal = N_thermal +Q_bars(:,:,k)*cte_mat_global(:,k)...
        *(z_top - z_bot);
    M_thermal = M_thermal +Q_bars(:,:,k)*cte_mat_global(:,k)...
        *(z_top^2 - z_bot^2);
end

N_thermal = N_thermal*delta_T;
M_thermal = M_thermal/2*delta_T;

thermal_loading = [N_thermal;...
                   M_thermal];