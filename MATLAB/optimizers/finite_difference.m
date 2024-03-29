function gradient = finite_difference(objective,x,h)
gradient = zeros(length(x),1);
for i = 1:length(x)
    x_step = zeros(length(x),1);
    x_step(i) = h;
    gradient(i) = (objective(x_step + x) - objective(x))/h;
end
end