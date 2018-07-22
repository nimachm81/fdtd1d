# fdtd1d

## A parallel one dimensional FDTD code for Maxwell equations. 

Finite difference time domain (FDTD) is one of the most efficient numerical
techniques for solving electromagnetic problems. It can handle a vast number
of practical electromagnetic problems such as electromagnetic interaction with
non-linear, dispersive, anisotropic and time varying media.

This code demonstrates how to write a concurrent FDTD using modern C++ features.

To compile and run:

```
$ cd build 
$ cmake .. 
$ make     
$ ./fdtd1d 
```

To set the number of threads, change the following line in main.cc:

``` C++
  int num_threads(1);
```

To record the fields and visualize the output, 

``` C++
  fdtd.SetTheWriteToFileFlag(false);
```

and then compile and run. To visualize the fields run the python script 
`plotOutput.py` from terminal:

```
$ python3 plotOutput.py
```








