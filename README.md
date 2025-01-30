# RLZ-RePair

## Description

This software computes the RePair encoding of a target sequence file via its RLZ parse representation and the corresponding reference file.

## Algorithmic Workflow

## Prerequisites

- [CMake](https://cmake.org/) 3.15 or higher.
- GCC 9+
- C++17-compatible compiler.
- OpenMP
- Python3

## Getting Started

### Building the Project

```
git clone https://github.com/rvarki/RLZ-RePair.git
cd RLZ-RePair
mkdir build
cd build
cmake ..
make -j
```
### Running the project

After building the project, an executable named rlz will be created in the build directory. Run it with:
```
./rlz-repair -r [reference file] -s [file to compress] [options] 
```

### Travis Example 

1. To run Travis' example, run the following command from the build directory.

```
./rlz-repair -r ../data/travis/travis_ref.txt -s ../data/travis/travis_seq.txt
```

### Rahul's Example

1. To run Rahul's example, run the following command from the build directory.

```
./rlz-repair -r ../data/rahul/rahul_ref.txt -s ../data/rahul/rahul_seq.txt
```






