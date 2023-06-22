import "magma_functions.m" : TraceFormulaGamma0AL, TraceFormulaGamma0ALNew;

procedure testPariVSMagma(N,k : New := false, OnlyPrimes := false, Prec := 1000)
    if (New and (not OnlyPrimes)) then
	printf("Error! Only supports newsubspace for prime Hecke operators!\n");
	assert false;
    end if;
    cmd := Sprintf("./src/traceALbatch_sta %o %o %o %o %o %o", N, N+1, k, Prec,
		   OnlyPrimes select 1 else 0, New select 1 else 0);
    System(cmd);
    fname := Sprintf("./data/traces_%o_%o.m", k, N);
    r := Read(fname);
    al_start := Index(r, "tracesAL");
    r := r[al_start..#r];
    start := Index(r, "[") + 1;
    fin := Index(r, "]") - 1;
    r := r[start..fin];
    traces_AL := [StringToInteger(x) : x in Split(r, ",")];
    trace_func := New select TraceFormulaGamma0ALNew else TraceFormulaGamma0AL;
    val_list := OnlyPrimes select PrimesUpTo(Prec) else [0..Prec];
    traces_magma := [trace_func(n, N, k) : n in val_list];
    assert traces_magma eq traces_AL[1..#val_list];
    return;
end procedure;

num_tests := 10;
max_Nk2 := 10^5;
max_level := 10000;
Ns := [1..max_level];
printf "testing pari vs magma implementation... (N;k) = ";
for i in [1..num_tests] do
    N := Random(Ns);
    max_weight := Floor(Sqrt(max_Nk2/N));
    ks := [2..max_weight by 2];
    k := Random(ks);
    printf "(%o;%o),", N, k;
    testPariVSMagma(N,k);
end for;

printf "testing pari vs magma implementation on new subspaces... (N;k) = ";
for i in [1..num_tests] do
    N := Random(Ns);
    max_weight := Floor(Sqrt(max_Nk2/N));
    ks := [2..max_weight by 2];
    k := Random(ks);
    printf "(%o;%o),", N, k;
    testPariVSMagma(N,k : New, OnlyPrimes);
end for;
printf "\n";

