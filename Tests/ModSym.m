procedure testPariVSMS(N,k,n)
    cmd := Sprintf("./src/traceALbatch_sta %o %o %o", N, N+1, k);
    System(cmd);
    fname := Sprintf("./data/traces_%o_%o.m", k, N);
    r := Read(fname);
    return_cmd := Sprintf("return traces_%o, tracesAL_%o;", N, N);
    traces_full, traces_AL := eval(r cat return_cmd);
    C := CuspidalSubspace(ModularSymbols(N,k,1));
    al := AtkinLehner(C,N);
    T := HeckeOperator(C,n);
    assert traces_AL[n+1] eq Trace(al*T);
    return;
end procedure;

num_tests := 10;
max_Nk2 := 20000;
max_level := 5000;
Ns := [1..max_level];
printf "testing pari vs modular symbols... (N;k;n) = ";
for i in [1..num_tests] do
    N := Random(Ns);
    max_weight := Floor(Sqrt(max_Nk2/N));
    ks := [2..max_weight by 2];
    k := Random(ks);
    n := Random([1..1000]);
    printf "(%o;%o;%o),", N, k, n;
    testPariVSMS(N,k,n);
end for;
printf "\n";
