# SVMGrad
SVMgrad is a compact library used to evaluate the decision function of a
Gaussian RBF Kernel Support Vector Machine, as well as the its first and
Second Derivative.

```
y = sign(Gamma(x))
Gamma   = \sum_{i=1}^{N_sv}\alpha_iy_ik(x,x_i) + b 
DGamma  = \sum_{i=1}^{N_sv}-1/2\sigma^2\alpha_iy_ik(x,x_i)(x-x_i)
DDGamma = ...
```
The evaluated SVM model can be learned with any toolbox: libSVM, SMVlight, EnsembleSVM, etc.. as long as the user creates a 
simplified struct model. In MATLAB one should create the following struct instance:
```
model.nClass: # of Classes (2 for binary)
model.nSV   : Total # of Support Vectors
model.b     : Offset for classification function
model.sigma : Gaussian RBF kernel Width
model.alphas: Values for the Lagrangian multipliers per SVs  [nSV]
model.y     : Labels corresponding to SVs                    [nSV]
model.SVs   : Set of Support Vectors                         [DxnSV]
```

##Matlab Usage
In the ```./matlab/``` folder you can find an example script together with sample 2d datasets and learnt models (through libSVM). Make sure to add all subfolders to your current directory and run the following script:
```
2d_svmgrad_example.m
```
This should generate the following plot, which shows the svm ```Gamma``` being evaluated on a mesh of points as well as its gradient ```DGamma```:

<p align="center">
<img src="https://github.com/nbfigueroa/SVMGrad/blob/master/img/2d-gamma.png" width="700">
</p>

##C++ Installation
Extract the SVMGrad folder to any location.
```
$ cd <SVM_grad_root_dir>
$ mkdir build
$ cd build
$ ccmake ..
$ make
$ sudo make install
```

##C++ Usage
bla bla bla


Note: Part of this code was derived from the [ASVM library](https://github.com/epfl-lasa/A-SVM) implemented by Dr. Ashwini Shukla.
