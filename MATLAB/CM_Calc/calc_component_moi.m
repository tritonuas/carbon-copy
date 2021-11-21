function i_moment = calc_component_moi(mass_vec,height_vec,width_vec)
    i_moment = 1/12*[mass_vec,mass_vec,mass_vec].*(height_vec.^2+width_vec.^2);
end