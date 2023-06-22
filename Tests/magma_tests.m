import "magma_functions.m" : testRelationNewSubspacesTrivial, testInverseRelationNewSubspacesTrivial,
       testRelationNewSubspaces, testInverseRelationNewSubspaces, get_trace, TraceFormulaGamma0ALNew;

num_tests := 10;
max_Nk2 := 20000;
max_level := 5000;
Ns := [1..max_level];
printf "testing magma formulas vs modular symbols... (N;k;p) = ";
for i in [1..num_tests] do
    N := Random(Ns);
    max_weight := Floor(Sqrt(max_Nk2/N));
    ks := [2..max_weight by 2];
    k := Random(ks);
    p := Random(PrimesUpTo(1000));
    printf "(%o;%o;%o),", N, k, p;
    testRelationNewSubspacesTrivial(N, k);
    testInverseRelationNewSubspacesTrivial(N, k);
    testRelationNewSubspaces(N, k, p);
    testInverseRelationNewSubspaces(N, k, p);
    assert TraceFormulaGamma0ALNew(p, N, k) eq get_trace(N, k, p : New);
end for;
printf "\n";
