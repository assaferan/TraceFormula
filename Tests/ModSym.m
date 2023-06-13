procedure testPariVSMS(N,n)
    cmd := Sprintf("./src/traceALbatch_sta %o %o", N, N+1);
    System(cmd);
    fname := Sprintf("./data/traces_%o.m", N);
    r := Read(fname);
    return_cmd := Sprintf("return traces_%o, tracesAL_%o;", N, N);
    traces_full, traces_AL := eval(r cat return_cmd);
    C := CuspidalSubspace(ModularSymbols(N,2,1));
    n := Random([1..1000]);
    al := AtkinLehner(C,N);
    T := HeckeOperator(C,n);
    assert traces_AL[n+1] eq Trace(al*T);
    return;
end procedure;

num_tests := 10;
max_level := 5000;
Ns := PrimesUpTo(max_level);
printf "testing pari vs modular symbols... (N;n) = ";
for i in [1..num_tests] do
    N := Random(Ns);
    n := Random([1..1000]);
    printf "(%o;%o),", N, n;
    testPariVSMS(N,n);
end for;
printf "\n";
