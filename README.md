# carbon-copy

This program is for optimizing the design of aircraft using simple optimization
techniques. The initial intended purpose for this program is to opimize a UAV 
design that can be manufactured by the undergraduate students of 
Triton Unmanned Aerial Systems (Triton UAS), a student project-team at UCSD, 
to use in the SUAS competition hosted by the AUVSI Seafarer Chapter.

## Installation
Everything is written in MATLAB, so MATLAB is required to run it.

Airfoil optimization (Windows):
Additionally, the XFoil executables and XFoil folder should be added at 
carbon-copy/MATLAB/Airfoil for the proper operating system. You can also figure
this out from their inclusion in the .gitignore file.
XFoil download and licensing can be found here: 
https://web.mit.edu/drela/Public/web/xfoil/

Airfoil optimization (Linux):
xfoil can be downloaded with sudo apt install xfoil

If there is an error when running the airfoil optimization, make sure the
correct lines are uncommented for your operating system in:
carbon-copy/MATLAB/Airfoil/ClCdData.m   (lines 17-20)
carbon-copy/MATLAB/Airfoil/modify_airfoil.m   (lines 34-38)


## Authors
Andrew Fletcher
Allyson Chen
Kevin Vo
Erik Lovekin
Tsz-Wai Kwok
Kaelan Tan
