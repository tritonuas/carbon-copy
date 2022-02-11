function [AR,wing_area] = geometric_outputs(wing_chord,wing_span)
wing_area = wing_chord*wing_span;
AR = wing_span^2/wing_area;
end