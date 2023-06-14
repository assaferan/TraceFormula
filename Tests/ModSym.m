procedure testPariVSMS(N,n)
    cmd := Sprintf("./src/traceALbatch_sta %o %o", N, N+1);
    System(cmd);
    fname := Sprintf("./data/traces_%o.m", N);
    r := Read(fname);
    al_start := Index(r, "tracesAL");
    r := r[al_start..#r];
    start := Index(r, "[") + 1;
    fin := Index(r, "]") - 1;
    r := r[start..fin];
    traces_AL := [StringToInteger(x) : x in Split(r, ",")];
    C := CuspidalSubspace(ModularSymbols(N,2,1));
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
