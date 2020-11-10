# carbon-copy

This program is for optimizing the design of aircraft using simple optimization
techniques. The initial intended purpose for this program is to opimize a UAV 
design that can be manufactured by the undergraduate students of 
Triton Unmanned Aerial Systems (Triton UAS), a student project-team at UCSD, 
to use in the SUAS competition hosted by the AUVSI Seafarer Chapter.

## Installation
Everything is written in MATLAB, so MATLAB is required to run it.
Additionally, the XFoil executables and XFoil folder should be added at 
carbon-copy/MATLAB/Airfoil for the proper operating system. You can also figure
this out from their inclusion in the .gitignore file.
XFoil download and licensing can be found here: 
https://web.mit.edu/drela/Public/web/xfoil/

## Authors
Andrew Fletcher
Tsz-Wai Kwok: Developing the wing weight considering and splitting the solution
    methods into their own functions.
Kaelan Tan: Improving the airfoil comparison code and writing a script in MATLAB 
    to run airfoil analysis on airfoil dat files using XFoil.