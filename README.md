# ReconstructionForOrthesis
Check these two papers, since those are the features being implemented here (actually the FPFH):
* [Point Feature Histograms](http://ezproxy.uniandes.edu.co:8080/login?url=http://ieeexplore.ieee.org/document/4650967/?part=1)
* [Fast Point Feature Histograms](http://ieeexplore.ieee.org.ezproxy.uniandes.edu.co:8080/document/5152473/?reload=true)

The idea here is to implement the FPFH so that that information can be given to a registration algorithm, namely the [Fast Global Registration](https://github.com/IntelVCL/FastGlobalRegistration) method developed by Vladlen Koltun, for it the perform the reconstruction (which is usually called registration).

If you want to work in this project you need to install the functions of FGR from the repository of the previous link (lots of instruction there, and I definitely recommend to try it on a Linux or OSX machine, NOT a Windows one) and download the MATLAB functions from this repository. 
