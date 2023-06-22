// SetDebugOnError(true);
// SetVerbose("ModularSymbols", 3);

import "magma_functions.m" : get_trace;

procedure testPariVSMS(N,k,n : New := false)
    is_prime := IsPrime(n);
    if (is_prime) then
	n_idx := #PrimesUpTo(n);
    else
	if (New) then
	    printf("Error! Only supports newsubspace for prime Hecke operators!\n");
	    assert false;
	end if;
	n_idx := n + 1;
    end if;
    cmd := Sprintf("./src/traceALbatch_sta %o %o %o %o %o %o", N, N+1, k, n, is_prime, New);
    System(cmd);
    fname := Sprintf("./data/traces_%o_%o.m", k, N);
    r := Read(fname);
    al_start := Index(r, "tracesAL");
    r := r[al_start..#r];
    start := Index(r, "[") + 1;
    fin := Index(r, "]") - 1;
    r := r[start..fin];
    traces_AL := [StringToInteger(x) : x in Split(r, ",")];
    // C := CuspidalSubspace(ModularSymbols(N,k,1));
    // al := AtkinLehner(C,N);
    // T := HeckeOperator(C,n);
    assert traces_AL[n_idx] eq get_trace(N,k,n : New := New); // Trace(al*T);
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
    n := Random([1..N]);
    printf "(%o;%o;%o),", N, k, n;
    testPariVSMS(N,k,n);
end for;

printf "testing trace on new subspace... (N;k;p) = ";
for i in [1..num_tests] do
    N := Random(Ns);
    max_weight := Floor(Sqrt(max_Nk2/N));
    ks := [2..max_weight by 2];
    k := Random(ks);
    p := Random(PrimesUpTo(N));
    printf "(%o;%o;%o),", N, k, p;
    testPariVSMS(N,k,p : New);
end for;

printf "testing trace on new subspace for bad primes... (N;k;p) = ";
for i in [1..num_tests] do
    N := Random(Ns);
    max_weight := Floor(Sqrt(max_Nk2/N));
    ks := [2..max_weight by 2];
    k := Random(ks);
    p := Random(PrimesDivisors(N));
    printf "(%o;%o;%o),", N, k, p;
    testPariVSMS(N,k,p : New);
end for;

printf "\n";
