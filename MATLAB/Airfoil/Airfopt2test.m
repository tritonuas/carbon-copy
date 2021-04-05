Cl = 1;
airfoil_dir_name = 'airfoil_database\';
for GuessAirfoil = [0123, 2222, 3217, 1712, 1010, 1315, 2212, 2214, 2217, 2010, 2012, 2014, 2016, 2018]
    Airfopt2(Cl, GuessAirfoil,airfoil_dir_name)
end
