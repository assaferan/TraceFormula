# TraceFormula
 Efficient computation of traceforms vis trace formula on Atkin-Lehner spaces
 
 Relies on the formula first described in [SZ], following the version appearing in [P].
 Using PARI/GP, as well as some of the code there based on Henri Cohen's algorithm for trace forms on the whole space and its implementation.
 Specifically using the original code for caching class number computations.
 
 At the moment the main branch only deals with weight 2 and prime level.
 The composite_level branch deals with any level and any weight.
 Both only compute traces on the Kohnen subspace (the Atkin-Lehner operator W_N where N is the level)
 
 [P] Popa, Alexandru A., On the trace formula for Hecke operators on congruence subgroups, II. Res. Math. Sci. 5 (2018), no. 1, Paper No. 3, 24 pp.
 [SZ] Skoruppa, Nils-Peter; Zagier, Don, Jacobi forms and a certain space of modular forms. Invent. Math. 94 (1988), no. 1, 113–146. 
