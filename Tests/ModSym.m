// SetDebugOnError(true);
// SetVerbose("ModularSymbols", 3);

function get_trace(N, k, n : New := false)
    C := CuspidalSubspace(ModularSymbols(N,k,1));
    if New then
	C := NewSubspace(C);
    end if;
    al := AtkinLehner(C,N);
    T := HeckeOperator(C,n);
    return Trace(T*al);
end function;

procedure testInverseRelationNewSubspacesTrivial(N, k)
    ms := [d : d in Divisors(N) | N mod d^2 eq 0];
    trace_1 := &+[Integers() | MoebiusMu(m) * get_trace(N div m^2, k, 1) : m in ms];
    assert trace_1 eq get_trace(N, k, 1 : New);
end procedure;

procedure testRelationNewSubspacesTrivial(N, k)
    N_primes := [d : d in Divisors(N) | IsSquare(N div d)];
    trace_1 := &+[Integers() | get_trace(N_prime, k, 1 : New) : N_prime in N_primes];
    assert trace_1 eq get_trace(N, k, 1);
end procedure;

function TrivialContribution(N, k, p)
    N_primes_2 := [d : d in Divisors(N) | IsSquare(N*p div d) and ((N div d) mod p^3 ne 0) and (d mod p ne 0)];
    trace_2 := &+[Integers() | get_trace(N_prime, k, 1 : New) : N_prime in N_primes_2];
    trace_3 := 0;
    if (N mod p eq 0) then
	N_primes_3 := [d : d in Divisors(N div p) | IsSquare(N div (d*p))];
	trace_3 := &+[Integers() | get_trace(N_prime, k, 1 : New) :  N_prime in N_primes_3];
    end if;
    return p^(k div 2) * trace_3  - p^(k div 2 - 1)*trace_2;
end function;

procedure testRelationNewSubspaces(N, k, p)
    N_primes_1 := [d : d in Divisors(N) | IsSquare(N div d) and ((N div d) mod p ne 0)];
    trace_1 := &+[ Integers() | get_trace(N_prime, k, p : New) : N_prime in N_primes_1];
    assert trace_1 + TrivialContribution(N, k, p) eq get_trace(N, k, p);
end procedure;

procedure testInverseRelationNewSubspaces(N, k, p)
    ms := [d : d in Divisors(N) | (N mod d^2 eq 0) and (d mod p ne 0)];
    // N_primes := [d : d in Divisors(N) | IsSquare(N div d) and ((N div d) mod p ne 0)];
    trace := &+[Integers() | MoebiusMu(m)*(get_trace(N div m^2, k, p) - TrivialContribution(N div m^2, k, p)) : m in ms];
    assert trace eq get_trace(N, k, p : New);
end procedure;

/*
primes := PrimesUpTo(100);
Nk2_bound := 500;
for k in [2..12 by 2] do
    max_N := Floor(Nk2_bound / k^2);
    for N in [1..max_N] do
	print k,N;
	testRelationNewSubspacesTrivial(N, k);
	testInverseRelationNewSubspacesTrivial(N, k);
	for p in primes do
	    testRelationNewSubspaces(N, k, p);
	    testInverseRelationNewSubspaces(N, k, p);
	    assert TraceFormulaGamma0ALNew(p, N, k) eq get_trace(N, k, p : New);
	end for;
    end for;
end for;
*/

procedure testPariVSMS(N,k,n)
    cmd := Sprintf("./src/traceALbatch_sta %o %o %o", N, N+1, k);
    System(cmd);
    fname := Sprintf("./data/traces_%o_%o.m", k, N);
    r := Read(fname);
    al_start := Index(r, "tracesAL");
    r := r[al_start..#r];
    start := Index(r, "[") + 1;
    fin := Index(r, "]") - 1;
    r := r[start..fin];
    traces_AL := [StringToInteger(x) : x in Split(r, ",")];
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
