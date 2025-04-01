# RLZ-RePair

## Description

RLZ-RePair efficiently compresses a target sequence file using a smaller reference file by first applying RLZ compression, then performing RePair compression on the RLZ parse. The algorithm aims to minimize the number of replacements that would occur if RePair were run directly on the sequence file by instead encoding the reference file as much as possible, since RLZ parse entries directly reference it.

## Algorithmic Workflow

The algorithm first applies **RLZ compression** to a target sequence file using a given reference file. This produces a .rlz parse file containing (ref position, length) pairs that losslessly represent the original sequence.

Then, the **RePair compression** step reads the RLZ parse and reference into memory and processes RLZ parse in a systematic manner to achieve the RePair compression of the sequence file. 

The characters represented by the RLZ parse are stored in two types of phrases:

- **Non-explicit phrases**: Store the start and end positions of the reference interval corresponding to its RLZ pair.

- **Explicit phrases**: Store the characters in an uncompressed format. 

At the start of the RePair section, all pairs in the RLZ parse are represented as non-explicit phrases (unless the user specifies otherwise).

If the most frequently occurring bi-gram is fully contained within the non-explicit phrases, the replacement only needs to occur in the reference, minimizing the total number of operations.

If the most frequent bi-gram spans the **phrase boundary** between two consecutive non-explicit phrases, the end of the first phrase and the start of the second are cleaved, and an explicit phrase containing the removed bi-gram is created between them. If the bi-gram spans a non-explicit phrase and an explicit phrase, only the start or end character of the non-explicit phrase is removed and added to the explicit phrase.

If the most frequent bi-gram occurs in the reference and either the start or end character of a non-explicit phrase partially overlaps with this position on the reference, this is called a **source boundary** condition. In this case, the start or end character of the non-explicit phrase must be removed and added to an adjacent explicit phrase.

**Time and memory bottlenecks in the software arise from managing explicit phrases.**

**For optimal performance, choose a reference that is highly similar to the sequence file so that the RLZ phrases are long.**

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

### Small Example Usage

In this section, we demonstrate how to run the tool using a small toy dataset and explore the different options available to the user. We highlight three (non-exclusive) ways to execute the tool from the build directory.

1. Default Settings

```
./rlz-repair -r ../data/travis/travis_ref.txt -s ../data/travis/travis_seq.txt
```

2. Verbose Option (-v [0,1,2])
```
./rlz-repair -r ../data/travis/travis_ref.txt -s ../data/travis/travis_seq.txt -v 2
```
| Verbosity Level | Output Details |
| -v 0 (default)  | Size of the RLZ parse and RePair compression. |
| -v 1 | Phrase size statistics, memory usage at start, max pair chosen each iteration, timing information. |
| -v 2 | Prints out the phrases each iteration and more. |

Each verbosity level includes all the information from the preceding levels.

> [!WARNING]
> Do not run with the -v 2 option for larger files as it will significantly slow the code down.

3. Experimental: Minimum Phrase Length Option (-m [num])

```
./rlz-repair -r ../data/travis/travis_ref.txt -s ../data/travis/travis_seq.txt -v 2 -m 2
```

Ensures non-explicit phrases are at least 2 characters long to reduce fragmentation. Smaller non-explicit phrases are converted to explicit phrases.

All threee commands will produce the following files in the in the data/travis directory: `travis_seq.txt.rlz`, `travis_seq.txt.C`, and `travis_seq.txt.R`. 

- `.rlz` file: contains the (pos,len) entries produced by the RLZ code. 
- `.C` file contains the compressed sequence file after being RePaired. 
- `.R` file contains the RePair rules to decompress the `.C` file back to the original sequence file.

4.  Decompression

```
./rlz-repair -d ../data/travis/travis_seq.txt
```
In the data/travis directory, the decompressed file will be labeled as `travis_seq.txt.out`.

5. Verify correctness
```
diff ../data/travis/travis_seq.txt ../data/travis/travis_seq.txt.out
```

There should be no output from this command which indicates that the files are identical.

### DNA Example

This is a larger example dataset used to verify the correctness of the software and assist in debugging. We will demonstrate its usage while applying some of the options explained in the previous section.

1. Run RLZ-RePair 

```
./rlz-repair -r ../data/dna/dna_ref.txt -s ../data/dna/dna_seq.txt -v1
```

2. Decompression

```
./rlz-repair -d ../data/dna/dna_seq.txt 
```

3. Verify correctness
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