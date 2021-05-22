function ClVals = getClVals(b,w,v,Cr,Ct)
    n = 11;
    ClVals = zeros(n,2);
    liftVals = getLiftVals(b,w);
    rho = 1.225;
    for i = 0:(n-1)
        y = liftVals(i+1,1);
        Ly = liftVals(i+1,2);
        Cy = Cr-(Cr-Ct)/(b/2)*y;
        ClVals(i+1,1) = y;
        ClVals(i+1,2) = Ly/((1/2)*rho*v^2*Cy);
    end
end