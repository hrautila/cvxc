
### CVX in C

The Cvxc package is rather direct porting of Python CVXOPT solvers into C. Major part of CvxC is libcxvc
library that includes conelp, cpl and cp solvers.

CvxC links by default with libARMAS linear algebra package. You can find libARMAS in github.com/hrautila/armas
repository. Matrix operations in the solvers are encapsulated in thin layer that is quite easy to replace with
some other linear algebra package.

## Compiling and installing

1. Download, compile and install libARMAS. Very short instructions for installing it under
   /homedir/.local directory: 

```
$ ./bootstrap.sh
$ mkdir build
$ cd build; ../configure --prefix=~/.local --enable-notypenames 
$ make CFLAGS="-O3"; make install
```

2. Download CvxC package

```
$ ./bootstrap.sh
$ mkdir build
$ cd build
$ CPPFLAGS=-I~/.local/include LDFLAGS=-L~/.local/lib ../configure --prefix=~/.local
$ make; make install
```

