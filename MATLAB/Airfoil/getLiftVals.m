function liftVals = getLiftVals(b,w)
   %{
    liftVals takes wingspan (b) & weight (w) to output n ellipse
    coordinates representing individual lift vectors at a distance
    y from origin.
   %}
    n = 11;
    liftVals = zeros(n,2);
    a = (8*((w/2))/0.9)/(pi*b);
    for i = 0:(n-1)
        y = (b/2)/10*i;
        liftVals(i+1,1) = y;
        liftVals(i+1,2) = sqrt(a^2*(1-(y^2)/(b/2)^2));
    end
end