import "magma_functions.m" : TraceFormulaGamma0AL;

procedure testPariVSMagma(N,k)
    cmd := Sprintf("./src/traceALbatch_sta %o %o %o", N, N+1, k);
    System(cmd);
    fname := Sprintf("./data/traces_%o_%o.m", k, N);
    r := Read(fname);
    return_cmd := Sprintf("return traces_%o, tracesAL_%o;", N, N);
    traces_full, traces_AL := eval(r cat return_cmd);
    traces_magma := [TraceFormulaGamma0AL(n, N, k) : n in [1..1000]];
    assert traces_magma eq traces_AL[2..1001];
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
printf "\n";

