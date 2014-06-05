denovo-qassembler
================
`denovo-qassembler` is a software package designed to reconstruct quaispecies variants from a set of next-generation
sequencing reads without using a reference sequence.

* [Download][]
    * [Older Versions][]
* [Building][]
    * [Dependencies][]
    * [Running cmake][]
    * [Disabling log4cxx][]
* [Running][]
    * [Examples][]
* [Contact][]

Download
--------
You can download the source code for the current release of `denovo-qassembler` below:

[`denovo qassembler 1.0.3`](denovo-qassembler-1.0.3.tar.gz)

### Older Versions

[`denovo qassembler 1.0.1`](denovo-qassembler-1.0.1.tar.gz)

[`denovo qassembler 1.0.2`](denovo-qassembler-1.0.2.tar.gz)


Recent Changes
--------------
### 1.0.3
* Fixed abundance estimation bug: edge weights were not reset prior to estimating abundance.
* Don't use `log1p` from boost to calculate relative abundance, use `log` from `cmath`.

### 1.0.2
* Fixed edge weight bug: edge weights are stored in an unsigned container, when decreasing edge weights make sure to check that the amount to decrease by is less than the current edge weight.
* Sort starting point vertices before constructing paths so that reconstructed paths are reproducible between executions.

Building
--------
### Dependencies
`denovo-qassembler` uses [`cmake`](http://www.cmake.org/) as a dependency detection and build automation tool. Binaries for
`cmake` are available for most operating systems. `denovo-qassembler` uses the [Boost C++ Libraries](http://www.boost.org/)
extensively for graph construction as well as pointer management. Finally, `denovo-qassembler` can optionally use
[`log4cxx`](http://logging.apache.org/log4cxx/) to display logging information with high granularity.

To install all dependencies in Ubuntu, you can use `apt-get`:

    sudo apt-get install g++ build-essential cmake libboost-all-dev liblog4cxx10-dev

To install all dependencies in Fedora, you can use `yum`:

    sudo yum install gcc-c++ cmake boost-devel log4cxx-devel zlib-devel

To install all dependencies in CentOS or Scientific Linux, you can use `yum`:

    sudo yum install gcc-c++ cmake boost-devel zlib-devel # CentOS and Scientific Linux do not ship log4cxx

### Running `cmake`
One of the primary features of `cmake` is out-of-source builds, meaning that all binary files will be built outside
of the source directory. Out-of-source builds make for a cleaner development environment because binary files are not
scattered throughout your source directory.

To build `denovo-qassembler`, extract the source tarball, change into the `build` directory, run `cmake`, then `make`:

    tar xf denovo-qassembler-1.0.3.tar.gz
    cd denovo-qassembler-1.0.3/build
    cmake ..
    make

`denovo-qassembler` has a comprehensive test suite which can be executed prior to installation to ensure correctness
of the build. To test `denovo-qassembler`, run:

    make test

`denovo-qassembler` can be installed to the default system directories (i.e., `/usr/local/bin` and `/usr/local/lib`) by
running:

    sudo make install

#### Disabling `log4cxx`
If `log4cxx` does not ship with your operating system (or you wish to disable it), you can still get some logging information
from `denovo-qassembler` by specifying log severity levels at compile time.

To disable `log4cxx`, pass `cmake` a compile time:

    cmake -DUSE_LOG4CXX=off ..

If you would like to enable some level of logging, you can run:

    cmake -DUSE_LOG4CXX=off -DLOGGING_SEVERITY=i

where `0 <= i <= 2` and larger values for `i` result in more log output. All logging written by `denovo-qassembler` when 
`log4cxx` is disabled is written to `stderr`.

Running
-------
`denovo-qassembler` has a variety of input options which change the way that variants are reconstructed and the way that
relative abundance is calculated.

| Option                          | Description                                                                                  | Default        | Type    | Required? |
|-------------------------------- | -------------------------------------------------------------------------------------------- | -------------- | --------| --------: |
| `--help`                        | Prints out all options and their description.                                                | disabled       | Boolean | No        |
| `--input-sequences` *f*         | fasta/fastq file containing sequencing reads.                                                | N/A            | String  | Yes       |
|                                 | Files compressed with `gzip` are allowed.                                                    |                |         |           |
| `--kmer-size` *i*               | The *k*-mer size used to construct the de Bruijn graph (*k* must be odd)                    | 31             | Integer | No        |
| `--pre-hash`                    | Hash all reads prior to constructing the graph.                                              | disabled       | Boolean | No        |
|                                 | May improve graph construction performance when used with `--aggressive-edge-removal`.       |                |         |           |
| `--aggressive-edge-removal` *i* | Remove edges from constructed graphs whose edge weight is below *i*.                         | N/A            | Integer | No        |
| `--print-graphs`                | Print the graph structures in DOT format, suitable for rendering with `graphviz`.            | disabled       | Boolean | No        |
| `--graph-dir` *d*               | Write DOT formatted graph files to the specified directory *d*.                              | `graphs/`      | String  | No        |
|                                 | Directory is created if necessary.                                                           |                |         |           |
| `--minimum-bases` *i*           | Remove disconnected sub-graphs with a single vertex referring to fewer                       | N/A            | Integer | No        |
|                                 | than *i* significant base pairs.                                                             |                |         |           |
| `--sequences`                   | Print the sequences corresponding to the paths generated by the path construction algorithm. | N/A            | Boolean | No        |
| `--sequence-dir` *d*            | Write sequences to the specified directory *d*. Directory is created if necessary.           | `sequences/`   | String  | No        |
| `--path-method` *m*             | Specify the method for constructing paths through the de Bruijn graph, one of                | `proportional` | String  | No        |
|                                 | `proportional`, `markov`, or `random`.                                                       |                |         |           |
| `--epsilon` *e*                 | Specify maximum allowable difference *e* between proportionally similar edge weights.        | 0.01           | Double  | No        | 
| `--minimum-length` *i*          | Do not report sequences that have a length less than *i*.                                    | N/A            | Integer | No        |
| `--abundance-method` *m*        | Specify the method for estimating relative abundance,                                        | `markov-chain` | String  | No        |
|                                 | one of `markov-chain` or `forward-algorithm`.                                                |                |         |           |
| `--log-config` *f*              | Specify a custom `log4cxx` configuration file.                                               | N/A            | String  | No        |

#### Examples
Given a compressed `fasta` file containing reads called `reads.fna.gz`, you could print out the set of de Bruijn graphs constructed using a *k*-mer size of 127:

    denovo-qassembler --input-sequences reads.fna.gz --kmer-size 127 --print-graphs --graphs-dir read-graphs

Another example is to reconstruct a set of variants using the Proportional method, estimate relative abundance with the `markov-chain` method and remove edges with edge weight below 2:

    denovo-qassembler --input-sequences reads.fna.gz --kmer-size 127 --aggressive-edge-removal 2 --sequences --sequence-dir read-sequences --abundance-method markov-chain

Contact
-------
For more information about `denovo-qassembler`, please contact:

* Franklin Bristow: <fbristow@cs.umanitoba.ca>
* Gary Van Domselaar: <gary.van.domselaar@phac-aspc.gc.ca>
* Michael Domaratzki: <mdomarat@cs.umanitoba.ca>
