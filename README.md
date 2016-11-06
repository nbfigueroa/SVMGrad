# SVMGrad
SVMGrad is a simple standalone library used to evaluate an SVM decision function (from a pre-learned SVM), as well as its first (gradient) and second order derivative, which can be used to compute the Jacobian and Hessian matrices of the function.

Note: This code is based on the [ASVM library](https://github.com/epfl-lasa/A-SVM) provided by Dr. Ashwini Shukla.

##Installation
####LINUX
Extract the SVMGrad folder to any location.
```
$ cd <SVM_grad_root_dir>
$ mkdir build
$ cd build
$ ccmake ..
$ make
$ sudo make install
```

##Usage
bla bla bla


##ThirdParty
###Libsvm
This package uses the function svmtrain from Libsvm. It contains a modified
subset of the original libsvm source which additionally returns the indices 
of the chosen support vectors within the model. Full version of Libsvm is 
available at 
[http://www.csie.ntu.edu.tw/~cjlin/libsvm](http://www.csie.ntu.edu.tw/~cjlin/libsvm)
Please read the COPYRIGHT file before using Libsvm.
