function central_moment = calc_central_moment(component_moi,weight_vec,position_vec,quarter_chord)
    distance = position_vec - quarter_chord ;
    moment_matrix = component_moi + weight_vec*(distance).^2;
    central_moment = sum(moment_matrix);
end