# RLZ-RePair

## Description

This software computes the RePair encoding of a target sequence file via its RLZ parse representation and the corresponding reference file. 

## Algorithmic Workflow

The main goal of the algorithm is to get a RePair encoding of the sequence file (represented through its RLZ parse) through doing the encoding with the reference file as much as possible. 

The RLZ portion of the software applies RLZ compression to a target sequence file using a provided reference file. This process generates a parse file (.rlz) that stores (reference start position, length) pairs, which losslessly represent the original sequence file. The RLZ implementation used in this software streams the sequence file and always tries to find the longest matching substring in the reference. Typically the reference should be much smaller than the sequence file to encode.

The RePair component then reads the RLZ parse and reference file into memory. The (pos, len) pairs in the RLZ parse represent a collection of "non-explicit" phrases. Each non-explicit phrase stores the start and end positions of its corresponding interval in the reference, which can be directly derived from the (pos, len) pair. 

In this algorithm, if the most frequently occuring bi-gram (in the sequence file) is completely contained within the intervals represented by non-explicit phrases, then the amount of work is minimal. In this case, replacements only need to be made in the reference which in theory should be much lower than the number of replacements in the sequence file. 

The bottleneck occurs when the most frequently ocurring bi-gram is not fully contained within all the non-explicit phrases. Either the bi-gram spans two consecutive phrases (phrase boundary condition) or the bi-gram on the reference only partially overlaps with the ends of non-explicit phrasess (source boundary condition). For both scenarios, in order for the non-explicit phrases to maintain compliance with the reference, the end(s) of the non-explicit phrase have to be sheared and an "explicit phrase" of the cleaved off end(s) has to be created. Per its name, an explicit phrase stores the characters it maintains in an uncompressed form. The majority of the time and memory resource usage of this algorithm in practice stems from having to create and maintain the non-explicit phrases for RePair.

For more detailed information on the algorithm and software. See the paper.

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
git submodule update --init --recursive
mkdir build
cd build
cmake ..
make -j
```
## Running the project

After building the project, an executable named rlz-repair will be created in the build directory. Run it with:
```
./rlz-repair -r [reference file] -s [file to compress] [options] 
```

### Small Example

In this section, we will showcase how to run the tool using a small toy dataset and some of the different options available to the user. We will highlight three (non-exclusive) ways to run the tool.

1a: To run RLZ-RePair with its default settings, run the following command from the build directory.  

```
./rlz-repair -r ../data/travis/travis_ref.txt -s ../data/travis/travis_seq.txt
```

This command will produce the following files in the data/travis directory: `travis_seq.txt.rlz`, `travis_seq.txt.C`, and `travis_seq.txt.R`. The `.rlz` file contains the (pos,len) entries produced by the RLZ code. The `.C` file contains the compressed sequence file after being RePaired. The `.R` file contains the RePair rules to decompress the `.C` file back to the original sequence file.

Running the default command (1a), only minimal information was printed to the screen. However, we do provide different levels of verbosity (-v {0,1,2}) that provide more insight about what the software is doing in the background. 

- v 0 (default): Provides information about the RLZ parse size and RePair compression size.
- v 1: Information v {0} provided + information about the size of the non-explicit phrases at the beginning, memory usage of different data structures at the beginning, the most frequent pair chosen at each iteration, and timing information about sections of the code (not exclusive).
- v 2 : Information v {0,1} provides + prints out the phrases, reference and more information for each part of the code every iteration.

> [!WARNING]
> Do not run with the -v 2 option for larger files as it will significantly slow the code down.

1b: Rerun the 1a command with the -v 2 option to see what is happening to the phrases each iteration

```
./rlz-repair -r ../data/travis/travis_ref.txt -s ../data/travis/travis_seq.txt -v 2
```

This command should have printed out how the phrases and reference look each iteration among other information. The file types produced are the same as the 1a command.

From 1b command, you should have seen that all the phrases at the start were non-explicit phrases as expected. However, some of these phrase were quite small (some only covering 1 character in the reference). Small phrases are detrimental since they are more likely to cause source and phrase boundaries conditions to be met which cause the tool to use a lot of resources. If changing the reference is not possible, we have included an experimental option (-m) which allows the user to set the minimum non-explicit phrase length at the start. If a non-explicit phrase does not meet this length, it will be converted to an explicit phrase at the start.

1c: Rerun the 1b command with the -m 2 option and observe at the start that the phrases less than 2 characters have been converted to explicit phrases.

```
./rlz-repair -r ../data/travis/travis_ref.txt -s ../data/travis/travis_seq.txt -v 2 -m 2
```

The file types produced are the same as the 1a command.

2:  The deompression command should be run as follows:

```
./rlz-repair -d ../data/travis/travis_seq.txt
```

This command will give some more information about the RePair compression. In the data/travis directory, the decompressed file will be labeled as `travis_seq.txt.out`.

3: Compare the sequence file to the decompressed file
```
diff ../data/travis/travis_seq.txt ../data/travis/travis_seq.txt.out
```

There should be no output from this command which indicates that the files are identical.

### DNA Example

This is a larger example dataset that was used to help determine the correctness of the software and to debug it. We will use some of the options explained in the previous section

1: Run RLZ-RePair from the build directory

```
./rlz-repair -r ../data/dna/dna_ref.txt -s ../data/dna/dna_seq.txt -v1
```

2: Decompress the file

```
./rlz-repair -d ../data/dna/dna_seq.txt 
```

3: Check the correctness
```
diff ../data/dna/dna_seq.txt ../data/dna/dna_seq.txt.out
```

## Usage
```
usage: rlz-repair [-h] [-r REF_FILE] [-s SEQ_FILE] [-m MIN] [-t THREADS] [-v {0,1,2}] [-d COMPRESS_PREFIX] [--log LOG_FILE]

 ____  _     _____     ____      ____       _      
|  _ \| |   |__  /    |  _ \ ___|  _ \ __ _(_)_ __ 
| |_) | |     / /_____| |_) / _ \ |_) / _` | | '__|
|  _ <| |___ / /|_____|  _ <  __/  __/ (_| | | |   
|_| \_\_____/____|    |_| \_\___|_|   \__,_|_|_|             
                ver 1.0.0

Compute the RLZ + RePair compression of a target sequence file 
using a reference file.

options:
  -h, --help            show this help message and exit
  --log LOG_FILE        Log filename.

Compression Options:
  -r REF_FILE, --ref REF_FILE
                        The reference file to be used for compression.
  -s SEQ_FILE, --seq SEQ_FILE
                        The sequence file to compress.
  -m MIN, --min MIN     The minimum phrase length threshold for an RLZ phrase to be stored as a non-explicit phrase.
  -t THREADS, --threads THREADS
                        Number of threads for RLZ section (default 1).
  -v {0,1,2}, --verbosity {0,1,2}
                        Set verbosity level (0 = none, 1 = basic, 2 = detailed).

Decompression Options:
  -d COMPRESS_PREFIX, --decompress COMPRESS_PREFIX
                        Prefix of the compressed files to decompress.
```