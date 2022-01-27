function [AR,chord] = geometric_outputs(wing_area,wing_span)
AR = wing_span^2/wing_area;
chord = wing_area/wing_span;
end