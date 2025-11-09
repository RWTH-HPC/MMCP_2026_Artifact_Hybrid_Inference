# Installation Guide # {#ugInstallationGuide}

[TOC]

# 1. GIT 
**Git** is a distributed version control system. Distributed means, that there is not only one single central repository containing and providing all revisions, but also ones cloned/copied by the users and existing on local devices for example. GitLab serves as the main central repository. That is, each developer has to update its repository to get changes from the central one, and has to push the changes from its repository to update the central one.
Each repository (cloned or not) consists of different branches. It is possible that branches exist only in one repository.  

## Prepare your GitLab account  

**SSH**

* Open the following link: https://git.rwth-aachen.de/
* If you have finished registration and got your TIM account:
    * Choose DFN-AAI Single Sign-On 
    * Choose RWTH Aachen and log in with your TIM account
* For people who do not have a TIM account yet:
    * Please prepare your **GitHub account** and ask a colleague to add you to the group.
    * You can log in by **GitHub**
* Setting your ssh keys: 
    * Upper right corner -> settings -> ssh keys, more details [here](https://docs.gitlab.com/ee/user/ssh.html) 

**E-Mail**  

* Click: Upper right corner -> settings -> Emails
* Check: Which E-Mail address is registered (your aia or student rwth)
* If you like, you can add a further address, such that later commits are allowed if signed by this address    

**m-AIA access**   

* Open the link: https://git.rwth-aachen.de/aia/MAIA/Solver  
    * If you have access you can continue with the following steps
    * If not: Contact one of your supervisors or one of maintainers

## Prepare your local Git
We can check our git settings by using:  
~~~
git config --list
~~~
Afterwards we may need to change our user name and e-mail address:  
~~~
git config --global user.name "John Doe"
git config --global user.email johndoe@example.com
~~~

This e-mail address will be used as default for our commits. Hence, it should agree with one e-mail specified before in your GitLab account.
We may want to set a default core editor such as vim or emacs:
~~~
git config --global core.editor vim
~~~

## Further information

http://ldap2.aia.rwth-aachen.de/mediawiki-1.22.1/index.php5/ZFS:GIT-Transition#Git_in_details
https://docs.gitlab.com/ee/tutorials/make_your_first_git_commit.html

# 2. Installation 
## Downloading m-AIA

There are several ways to access the m-AIA source code. 
If you intend to be an developer, you are recommended to clone the repository.

**Cloning**

* Hosts except Hawk:
you can just use:
~~~
git clone git@git.rwth-aachen.de:aia/MAIA/Solver.git ./ 
~~~
Tips: this is the way only if you have set your ssh key on gitlab. You can also clone with HTTPS that will use your username and password of GitLab.

* Host Hawk:
Even if you have access to GitLab, the HAWK systems can only connect to certain approved IPs, you can't use git like you are used to on AIA:
On the HAWK systems git has to connect to the git-servers of the RWTH through an ssh-tunnel or HTTP tunnel.
For more details, you can check [here](http://ldap2.aia.rwth-aachen.de/mediawiki-1.22.1/index.php5/Com:hpc#Git_access)

Following successful cloning, an individual copy of the whole remote repository has been created on your local device.

## Building
### Requirements

**Python3.x and CMAke**

m-AIA requires Python 3 to compile and run the Python scripts. 
Make sure you have the python3 executables set up correctly 
in your environment.CMAke is used to control the compilation process, and version restriction is **CMAke >= 2.8**.

**Compilers**

The following compilers with version are currently supported and tested:
* GCC >= 9.1
* Clang >= 9
* Intel >= 19.1
* PGI >= 19.2
* NVHPC >=21.11

**Packages**

Those software packages are required for installation:

* Parallel NetCDF >= 1.9
* OpenMPI >= 4.0
* FFTW (with MPI support) >= 3.3.2
* HDF5

The automatic configuration files are provided for some hosts, including:
* AIA
* JUWELS (Forschungszentrum Juelich)
* JURECA (Forschungszentrum Juelich)
* HAWK (HLRS Stuttgart)
* Claix (RWTH Aachen)
* Klogin
* Mac
* PizDaint
* Power8
* PSG

First the `configure.py` script prepares the compilation process for your specific host and desired
compiler settings. This assumes that the necessary information for your system is available, which
is the case for the afore-mentioned hosts. 
New hosts are allowed to be added following the instructions here:`~/Solver/auxiliary/hosts/README`. You can also find a file named 
*Host.cmake.in* whick is a reference for new hosts' <em><HOST_NAME>.cmake</em>. After getting all configuration files, the configuration can be continued.

To check for all available compilers and build types on your current host, run `./configure.py`.
You can then choose one by running, e.g.
~~~
./configure.py 1 1
~~~
or
~~~
./configure.py GNU debug
~~~
which will set up the Makefile for compilation with GNU compiler in debug mode. Furthermore,
`./configure.py -h` can be used to check for other optional arguments which will be helpful if you have special needs.
For instance, `./configure.py 1 1 --without-hdf5` indicates this configuration does not link against HDF5 library and disables 
HDF5-related code.

Calling the configure.py script with the arguments `COMPILER` and `BUILDTYPE` will then create a new directory of the format `build_COMPILER_BUILDTYPE/`, which contains all the files needed for building MAIA. Compilation is initiated via the Makefile in the MAIA root directory, by calling
~~~
make <build_COMPILER_BUILDTYPE>
~~~
and supplying optional arguments like `-j NO_NODES`. Since a separate directory is created for each build configuration, multiple instances of MAIA (say GNU/debug and GNU/production) can be maintained at the same time. Calling make with no targets (which defaults to `all`) or with a valid target, that is not a directory (like `all`, `coverage`, `depend`, `what` etc...), each applicable build directory will be treated with the target supplied via the command line. Individual configurations, e.g. `build_gnu_production`, can be targeted via `make build_gnu_production/TARGET` (with `TARGET` being a valid target for the build configuration, like `all`, `clean`, `coverage`, `what`).
Additionally, a new target `link` is available, which calls for each previously configured build config `make` with the target `all` and then symlinks the executables of each build config to the maia root directory in the form `maia_COMPILER_BUILDTYPE`. If it is called directly on an individual build directory via e.g. `make build_gnu_production/link`, it only creates a link to that executable in the MAIA root dir and names it `maia` and keeps the original name. This way one can maintain one default or main version of the solver.
Even more customization is possible by the new optional argument `--build-dir-name BUILDDIRNAME` in the configure.py script. This way, compilation settings like `--with-cantera` or `--without-hdf5` can be distinguished even if compiler and build type chosen are the same.
For cleanup, use `make distclean` to restore the pre-configuration state.


**Libraries installation**

On Hosts outside the AIA, you will have to provide additional external libraries. You can use `module spider fftw` to check if this
module already exists. Then you can `module load fftw` to load the library. If there is not, then you have to install it manually.
It is recommended to locate your libraries in the directory $HOME/libraries/ and use folder names without version numbers ( like $HOME/libraries/netcdf/ ). 
You can also use a link pointing to that directory, if you like, for example, if you have installed different versions of an external library.
For example, to install the libraries FFTW, Parallel-Netcdf and Netcdf
~~~
tar -zxvf nameOfTheLib.tar.gz
cd nameOfTheLib
./configure --prefix=PathToTheLib(use_pwd) CC=nameOfC-Compiler CXX=nameOfC++-Compiler F77=nameOfF77-Compiler \
   F90=nameOfF90-Compiler MPICC=nameOfMPI-C-Compiler MPICXX=nameOfMPI-C++-Compiler MPIF77=nameOfMPI-F77-Compiler \ 
   MPIF90=nameOfMPI-F90-Compiler
make
make check
make install
~~~

After configuring successfully, then to start compiling, run: 
~~~
make -j <np>
~~~
where <np> is the number of cores used to compile the code. 

when you get executable file *maia* in folder called *src* , you can use `maia -h` to get a full list of command line options like:
~~~
  -b, --build-type         show build type that was used
  -c, --compiler           show used compiler
  -d, --debug LEVEL        use specified debug level
  -h, --help               this help screen
  -i, --input DIRECTORY    specified directory is used for all input files
  -o, --output DIRECTORY   specified directory is used for all output files
  -v, --version            show MAIA version
~~~
This indicates that the compilation has ended successfully, m-AIA is ready to use.
