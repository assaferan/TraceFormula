# TraceFormula
 Efficient computation of trace formula on Kohnen's plus-space of classical modular forms S_k(N)^+
 
 ## Requirements
 Pari/GP
 
 ## Installation
 ```
 > ./autogen.sh
 > ./configure
 > ./make
 ```
 
 ## Running
 ```
 > ./src/traceALbatch_dyn <from> <to> <k>
 ```
 
 Creates files ```./data/traces_<k>_<N>.m``` for N in the range ```[<from>, <to>)```.
 Each files contains a two lists of traces - one of the operators T_n and one of the operators T_n * W_N for n up to the Sturm bound, or up to 1000, if the Sturm bound is lower, on the space S_k(N).
