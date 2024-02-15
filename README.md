# TraceFormula
 Efficient computation of traceforms vis trace formula on Atkin-Lehner spaces
 
 Relies on the formula first described in [SZ], following the version appearing in [P].
 Using PARI/GP, as well as some of the code there based on Henri Cohen's algorithm for trace forms on the whole space and its implementation.
 Specifically using the original code for caching class number computations.

 ## Requirements

- [autotools](https://www.gnu.org/software/automake/manual/html_node/Autotools-Introduction.html)
- [gcc](https://gcc.gnu.org/)

This code has only been tested so far on:

macOS Monterey 12.3.1, compiled with clang 12.0.0 on Intel i9 8-Core Processor.

macOS Catalina 10.15.7, compiled with clang 12.0.0 on Intel Core i5 Quad-Core Processor.

linux Ubuntu 22.04.2, compiled with gcc 11.3.0 on AMD Ryzen ThreadRipper 2970WX 24-Core Processor.

The main (and develop) branch require the following libraries:

- [libgmp-dev](https://gmplib.org)
- [libmpfr-dev](https://www.mpfr.org)
- [libtool-bin](https://www.gnu.org/software/libtool/)
- [libpari-dev](https://pari.math.u-bordeaux.fr)

Installation of all the requirements on either linux or macOS can be done as follows.

    brew update && brew install autoconf automake libtool libgmp-dev libmpfr-dev libtool-bin libpari-dev

or
    
    sudo apt-get update && sudo apt-get install autoconf automake libtool libgmp-dev libmpfr-dev libtool-bin libpari-dev

## Installation

To build as a C library on a standard Linux/Mac OS system with autotools:

    ./autogen.sh
    ./configure --prefix=DIR
    make
    make install

Or all at once:

   ./autogen.sh && ./configure && make && make install

## Usage

The executable is src/traceALbatch_dyn.
There is also a static version of the code named src/traceALbatch_sta.
It runs with command-line arguments as follows

src/traceAL_dyn <from> <to> <weight> (<num_traces>) (<only_primes>) (<new>)

where the  arguments are:

<from> - The starting level for which to compute traces.

<to> - The final level for which to compute traces.

<weight> k - the weight of modular forms to compute the trace for.

<num_traces> The number of traces to compute. If not specified computes up to max(B, 30sqrt(p), 1000) where B is the Sturm bound.

<only_primes> Either 0 or 1. If specified as 1, only computes traces for Hecke operators at primes T_p.

<new> Either 0 or 1. If specified as 1 computes traces on spaces of newforms. 

Output is written to files in the subfolder data with names traces_k_N.m where k is the weight and N is the level
    
Example run:

src/traceALbatch_sta 1 100 2
 
 [P] Popa, Alexandru A., On the trace formula for Hecke operators on congruence subgroups, II. Res. Math. Sci. 5 (2018), no. 1, Paper No. 3, 24 pp.
 
 [SZ] Skoruppa, Nils-Peter; Zagier, Don, Jacobi forms and a certain space of modular forms. Invent. Math. 94 (1988), no. 1, 113–146. 
