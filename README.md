# Welcome to the SWASH git repository

SWASH is a general-purpose numerical tool for simulating unsteady, non-hydrostatic, free-surface, rotational flow and transport phenomena in coastal waters as driven by waves,
tides, buoyancy and wind forces. It provides a general framework for describing wave transformations from deep water to a beach, port or harbour, complex changes to rapidly varied
flows, and density driven flows in coastal seas, estuaries, lakes and rivers.

For more in-depth background and scientific documentation, the reader is referred to the [SWASH website](https://swash.sourceforge.io).
Please also check the [release notes](https://swash.sourceforge.io/modifications/modifications.htm) for any additional information on the current version **11.01**.

This Readme provides a brief overview of software installation and configuration instructions for users and developers.
Please see the [Implementation Manual](https://swash.sourceforge.io/online_doc/swashimp/swashimp.html) for additional documentation.

In addition to the installation, a brief outline on how to run the model is given below.

The SWASH software can be used freely under the terms of the [GNU General Public License](https://gitlab.tudelft.nl/citg/wavemodels/swash/-/blob/main/LICENSE).
It is permitted to copy, reuse, adapt and distribute the SWASH source code provided that proper reference is made to the original work.

## installation

#### prerequisites

To install SWASH on your local system, CMake and Ninja (or GNU make) need to be installed first.
We recommend to use [CMake 3.12+](https://cmake.org/) for building SWASH (check the version by typing `cmake --version`).
CMake is a build system and makes use of scripts (or configuration files) that control the build process.
There are installers available for Windows, Linux and macOS. See the
[download](https://cmake.org/download/) page for CMake installation instructions.

[Ninja](https://ninja-build.org/) is one of the many build tools to create executable files and libraries from source code.
The way it works is very similar to GNU make (or NMAKE for Windows); for example, it does not rebuild things that are already up to date.
Ninja can be downloaded from its [git repository](https://github.com/ninja-build/ninja/releases).

In addition to the build tools, a Perl package must be available on your local computer.
Usually, it is available for macOS, Linux and a UNIX-like operating system. Check it by typing `perl -v`.
You can download Perl for MS Windows from [Strawberry Perl](https://strawberryperl.com).
The Perl version should be at least 5.0.0 or higher.

Finally, SWASH also requires a Fortran90 compiler to be present in your environment.
Popular Fortran compilers are [gfortran](https://gcc.gnu.org/fortran/) and
[Intel<sup>&reg;</sup> Fortran Compiler Classic](https://www.intel.com/content/www/us/en/developer/articles/tool/oneapi-standalone-components.html#fortran)
(as part of the Intel<sup>&reg;</sup> oneAPI HPC Toolkit) and both support the OpenMP standard.
Please check this [page](https://fortran-lang.org/learn/os_setup/install_gfortran) for the installation of gfortran on your platform.

Due to the use of ANSI standard Fortran90 the SWASH source code can be ported to various architectures (e.g., Windows, Linux, macOS and Unix-like systems).
Currently, the build scripts support the following Fortran compilers:

1. GNU
1. Intel<sup>&reg;</sup>
1. Portland Group
1. Lahey
1. IBM XL Fortran

#### instructions

##### 1. clone the repo and navigate to the top level source directory

```bash
$ git clone https://gitlab.tudelft.nl/citg/wavemodels/swash.git && cd swash
```

##### 2. create the build directory

At the top of SWASH source directory execute the following commands

```bash
$ mkdir build && cd build
```

This step is required to perform an out-of-source build with CMake, that is, build files will not be created in the `/swash/src` directory.

##### 3. build the software

Two CMake configuration files are provided as required for the build. They are placed in the following source directories: `./swash/CMakeLists.txt` and `./swash/src/CMakeLists.txt`.

The following two CMake commands should suffice to build SWASH

```bash
$ cmake .. -G Ninja
$ cmake --build .
```

The first command refers to the source directory where the main configuration file is invoked. The second command carries out the building in the build directory.

The package is actually built by invoking Ninja. An alternative would be to use GNU make, as follows

```bash
$ cmake .. -G "Unix Makefiles"
$ make
```

or just (in case your OS is Unix-like)

```bash
$ cmake ..
$ make
```

However, we recommend Ninja because it is faster than GNU make.

##### 4. install the package

To install SWASH, run either

```bash
$ cmake --install .
```

or with the GNU make

```bash
$ make install
```

The default install directory is `/usr/local/swash` (Unix-like operating systems, including macOS) or `C:\Program Files\swash` (Windows).
Instead, you may install SWASH in any other user-defined directory, as follows

```bash
$ cmake --install . --prefix /somewhere/else/other/than/default/directory
```

or

```bash
$ cmake .. -DCMAKE_INSTALL_PREFIX=/somewhere/else/other/than/default/directory
$ make install
```

After installation a number of subdirectories are created.
The executables end up in the `/bin` directory, the archive/library files in `/lib`, and the module files in `/mod`.
Additionally, the `/doc` folder contains the pdf documents, the folder `/tools` consists of some useful scripts and the `/misc` directory
contains all of the files that do not fit in other folders (e.g., a machinefile and an edit file).

Please note that the installation can be skipped (though not recommended). Executables and libraries are then located in subdirectories of the build directory.

#### configuring the build

The build can be (re)configured by passing one or more options to the CMake command with prefix `-D`. A typical command line looks like

```bash
$ cmake .. -D<option>=<value>
```

where `<value>` is a string or a boolean, depending on the specified option. The table below provides an overview of the non-required options that can be used.

|  option                  | value type |               description                 | default value           |
|:------------------------:|:-----------|:------------------------------------------|:-----------------------:|
| `CMAKE_INSTALL_PREFIX`   | string     | user-defined installation path            | `/usr/local/swash`      |
| `CMAKE_Fortran_COMPILER` | string     | full path to the Fortran compiler         | determined by CMake     |
| `MPI`                    | boolean    | enable build with MPI                     | `OFF`                   |
| `CMAKE_VERBOSE_MAKEFILE` | boolean    | provide verbose output of the build       | `OFF`                   |

For example, the following commands

```bash
$ cmake .. -GNinja -DMPI=ON
$ cmake --build .
```

will configure SWASH to be built created by Ninja that supports parallel computing using the MPI paradigm.
Note that CMake will check the availability of MPI libraries within your environment.

The system default Fortran compiler (e.g., f77, g95) can be overwritten as follows

```bash
$ cmake .. [options] -DCMAKE_Fortran_COMPILER=/path/to/the/desired/compiler/including/the/name/of/compiler
```

Finally, if CMake fails to configure your project, then execute

```bash
$ cmake .. [options] -DCMAKE_VERBOSE_MAKEFILE=ON
```

which will generate detailed information that may provide some indications to debug the build process.

#### clean up the build files

To remove the build directory and all files that have been created after running `cmake --build .`, run at the top level of your project the following command:

```bash
$ cmake -P clobber.cmake
```

(The `-P` argument passed to CMake will execute a script *\<filename\>.cmake*.)

## getting started

*Note: before start using the SWASH package, it is suggested to first read
Chapters [2](https://swash.sourceforge.io/online_doc/swashuse/node2.html#ch:defin) and
[3](https://swash.sourceforge.io/online_doc/swashuse/node21.html#ch:inout)
of the [SWASH's User Guide](https://swash.sourceforge.io/online_doc/swashuse/swashuse.html).
Also, Chapter [5](https://swash.sourceforge.io/online_doc/swashuse/node41.html#ch:yourcom) is recommended.*

#### run modes

There are two different modes in which you can run SWASH:

1. **serial**: for laptops using one processor (recommended for small tests)
1. **parallel, distributed (MPI)**:  for distributed memory systems (laptop/desktop using multiple processors, clusters, HPC)

See the above installation instructions for getting the proper executable.

#### how to run

The general run procedure is as follows:

1. complete or modify your command file `INPUT`
1. run the SWASH model:
   ```bash
   $ ./swash.exe
   ```
1. check the created `PRINT` file for warning and error messages
1. repeat if needed

For faster simulation on a cluster, replace the run command by

```bash
$ mpirun -np <n> swash.exe
```

with `<n>` the number of desired nodes.

*Note: in the current version (v 11.01), the use of unstructured mesh is not supported in parallel mode.*

The above procedure can be done automatically using the script `/bin/swashrun` (or `\bin\swashrun.bat` in case of Windows), provided that the
environment variable `PATH` has been adapted by including the path of the `/bin` directory.

For more details, consult the [Implementation manual](https://swash.sourceforge.io/online_doc/swashimp/node12.html).

## documents

See
1. the [Implementation Manual](https://swash.sourceforge.io/online_doc/swashimp/swashimp.html) that describes in detail the installation and the usage of the SWASH model
1. the [User Manual](https://swash.sourceforge.io/online_doc/swashuse/swashuse.html) that provides the specifications for the input of the SWASH model
1. the [Scientific/technical documentation](https://swash.sourceforge.io/online_doc/swashtech/swashtech.html) that discusses the physical and mathematical details and the discretizations that have been implemented in the SWASH software
1. a bunch of scientific documents that can be found [here](https://swash.sourceforge.io/online_doc/online_doc.htm)

## bugs and questions

For bug reports please send to the [SourceForge mailing list](http://sourceforge.net/mail/?group_id=385225).


<small>&copy; Copyright 2025  Marcel Zijlema</small>
