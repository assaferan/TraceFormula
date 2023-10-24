function Sfast(N, u, t, n)
  fac := Factorization(N*u);
  num_sols := 1;
  for f in fac do
    p,e := Explode(f);
    if p eq 2 then
      e2 := Valuation(1-t+n, 2);
      if (e2 eq 0) or (e gt e2) then
        return 0;
      end if;
      if IsEven(t) then
        num_sols *:= 2^(e-1);
      end if;
    else
      y := t^2-4*n;
      e_y := Valuation(y, p);
      if e_y lt e then
        if IsOdd(e_y) then
          return 0;
        end if;
	y_0 := y div p^e_y;
        is_sq := IsSquare(Integers(p)!y_0);
        if not is_sq then
	  return 0;
        end if;
        num_sols *:= p^(e_y div 2) * 2;
      else
	num_sols *:= p^(e div 2);  
      end if;
    end if;
  end for;
  return num_sols div u;
end function;

function S(N, u, t, n)
  assert N mod u eq 0;
  assert (t^2 - 4 * n) mod u^2 eq 0;
  return [x : x in [0..N-1] | GCD(x,N) eq 1 and (x^2 - t*x + n) mod (N*u) eq 0];
end function;

function phi1(N)
  primes := [f[1] : f in Factorization(N)];
  if IsEmpty(primes) then return N; end if;
  return Integers()!(N * &*[1 + 1/p : p in primes]);
end function;

function B(N, u, t, n)
// assert Sfast(N,u,t,n) eq #S(N,u,t,n);
    return #S(N,u,t,n) * phi1(N) div phi1(N div u);
  // return Sfast(N,u,t,n) * phi1(N) div phi1(N div u);  
end function;

function C(N, M)
  a, b, c, d := Explode(M);
  G_M := GCD(c, d-a, b);
  return B(N, GCD(G_M, N), Trace(M), Determinant(M));
end function;

function C(N, u, t, n)
  return &+[B(N, u div d, t, n) * MoebiusMu(d) : d in Divisors(u)];
end function;

function Lemma4_5(N, u, D)
    assert N mod u eq 0;
    assert D mod (u^2) eq 0;
    ret := 1;
    fac := Factorization(N);
    for pa in fac do
	p, a := Explode(pa);
	i := Valuation(u, p);
	if (i eq 0) then continue; end if;
	b := Valuation(D, p);
	if IsOdd(p) then
	    if (i eq a) then ret *:= p^(Ceiling(a/2)); continue; end if;
	    if ((i le b-a) and IsEven(a-i)) then
		ret *:= (p^Ceiling(i/2) - p^(Ceiling(i/2)-1));
	    elif ((i eq b-a+1) and IsEven(a-i)) then
		ret *:= - p^(Ceiling(i/2)-1);
	    elif ((i eq b-a+1) and IsOdd(a-i)) then
		ret *:= p^Floor(i/2) * LegendreSymbol(D div p^b, p);
	    else
		return 0;
	    end if;
	else // p = 2
	    if (i eq a) then
		if (b ge 2*a+2) or ((b eq 2*a) and (((D div 2^b) mod 4) eq 1)) then 
		    ret *:= p^(Ceiling(a/2));
		elif (b eq 2*a+1) or ((b eq 2*a) and (((D div 2^b) mod 4) eq 3)) then
		    ret *:= -p^(Ceiling(a/2)-1);
		end if;
		continue;
	    end if;
	    if ((i le b-a-2) and IsEven(a-i)) then
		// print "case 1";
		ret *:= p^(Ceiling(i/2)-1);
	    elif ((i eq b-a-1) and IsEven(a-i)) then
		// print "case 2";
		ret *:= - p^(Ceiling(i/2)-1);
	    elif ((i eq b-a) and IsEven(a-i)) then
		// print "case 3";
		ret *:= p^(Ceiling(i/2)-1) * KroneckerCharacter(-4)(D div p^b);
	    elif ((i eq b-a+1) and IsOdd(a-i) and ((D div p^b) mod 4 eq 1) ) then
		// print "case 4";
		ret *:= p^Floor(i/2) * KroneckerSymbol(D div p^b, p);
	    else
		// print "returning 0";
		return 0;
	    end if;
	end if;
    end for;
    return ret;
end function;

function Cfast(N, u, t, n)
    S := [x : x in [0..N-1] | (GCD(x,N) eq 1) and (((x^2 - t*x + n) mod N) eq 0)];
    return #S * Lemma4_5(N, u, t^2 - 4*n);
end function;

function Hurwitz(n)
    if n eq 0 then
	return -1/12;
    end if;
    t_sum := 0;
    
    for d in Divisors(n) do
	is_sq, f := IsSquare(d);
	if is_sq then
            D := -n div f^2;
	    if D mod 4 in [0,1] then
		O := QuadraticOrder(BinaryQuadraticForms(D));
		h := #PicardGroup(O);
		w := #TorsionSubgroup(UnitGroup(O));
		t_sum +:= (h/w);
	    end if;
	end if;
    end for;
    return 2*t_sum;
end function;

function H(n)
  if n lt 0 then
    is_sq, u := IsSquare(-n);
    return (is_sq select -u/2 else 0);
  end if;
  if n eq 0 then
    return -1/12;
  end if;
  if n mod 4 in [1,2] then
    return 0;
  end if;

  ret := &+[ClassNumber(-n div d) : d in Divisors(n)
		       | IsSquare(d) and (n div d) mod 4 in [0,3] ];
  if IsSquare(n) and IsEven(n) then
    ret -:= 1/2;
  end if;
  if n mod 3 eq 0 and IsSquare(n div 3) then
    ret -:= 2/3;
  end if;
  return ret;
end function;

function Phi(N, a, d)
  ret := 0;
  for r in Divisors(N) do
    s := N div r;
    // scalar := r eq s select 1/2 else 1;
    g := GCD(r,s);
    if GCD(N, a-d) mod g eq 0 then
      alpha := CRT([a,d],[r,s]);
      if (GCD(alpha,N) eq 1) then
        ret +:= EulerPhi(g);
      end if;
    end if;
  end for;
  return ret;
end function;

// Gegenbauer polynomials
function P(k, t, m)
  R<x> := PowerSeriesRing(Rationals(), k-1);
  return Coefficient((1 - t*x+m*x^2)^(-1), k-2);
end function;

// This should yield -2*(A1 + A2)
function S1Popa(n,N,k)
    S1 := 0;
    max_abst := Floor(SquareRoot(4*n));
    for t in [-max_abst..max_abst] do
	for u in Divisors(N) do
	    if ((4*n-t^2) mod u^2 eq 0) then
		// print "u = ", u, "t = ", t;
		S1 +:= P(k,t,n)*H((4*n-t^2) div u^2)*C(N,u,t,n);
		// assert H((4*n-t^2) div u^2) eq Hurwitz((4*n-t^2) div u^2);
		// S1 +:= P(k,t,n)*Hurwitz((4*n-t^2) div u^2)*C(N,u,t,n);
		// print "S1 = ", S1;
	    end if;
	end for;
    end for;
    return S1;
end function;

// This should yield -2*A3
function S2Popa(n,N,k)
    S2 := 0;
    for d in Divisors(n) do
	a := n div d;
	S2 +:= Minimum(a,d)^(k-1)*Phi(N,a,d);
    end for;
    return S2;
end function;

function TraceFormulaGamma0(n, N, k)
    S1 := S1Popa(n,N,k);
    S2 := S2Popa(n,N,k);
    ret := -S1 / 2 - S2 / 2;
    if k eq 2 then
	ret +:= &+[n div d : d in Divisors(n) | GCD(d,N) eq 1];
    end if;
    return ret;
end function;

function PhiAL(N, a, d)
  return EulerPhi(N) / N;
end function;

function Phil(N, l, a, d)
    l_prime := N div l;
    ret := 0;
    for r in Divisors(l_prime) do
	s := l_prime div r;
	g := GCD(r,s);
	if (((a-d) mod g eq 0) and (GCD(a,r) eq 1) and (GCD(d,s) eq 1)) then
	    ret +:= EulerPhi(g);
	end if;
    end for;
    return EulerPhi(l) * ret / l;
end function;

function alpha(n)
    fac := Factorization(n);
    ret := 1;
    for fa in fac do
	if (fa[2] in [1,2]) then ret := -ret; end if;
	if (fa[2] ge 4) then return 0; end if;
    end for;
    return ret;
end function;

// Now this seems to work - tested for N <= 100, even k, 2<=k<=12 and 1 <= n <= 10, n = N
// This formula follows Popa - On the Trace Formula for Hecke Operators on Congruence Subgroups, II
// Theorem 4. 
// (Also appears in Skoruppa-Zagier, but this way of stating the formula was easier to work with).
function TraceFormulaGamma0AL(n, N, k)
    if (n eq 0) then return 0; end if; // for compatibility with q-expansions
  S1 := 0;
//  max_abst := Floor(SquareRoot(4*N*n));
  max_abst := Floor(SquareRoot(4*N*n)) div N;
  // ts := [t : t in [-max_abst..max_abst] | t mod N eq 0];
  // for t in ts do
  for tN in [-max_abst..max_abst] do
    t := tN*N;
    for u in Divisors(N) do
      if ((4*n*N-t^2) mod u^2 eq 0) then
	S1 +:= P(k,t,N*n)*H((4*N*n-t^2) div u^2)*C(1,1,t,N*n)
	       *MoebiusMu(u) / N^(k div 2-1);
      end if;
    end for;
  end for;
  S2 := 0;
  for d in Divisors(n*N) do
    a := n*N div d;
    if (a+d) mod N eq 0 then 
      S2 +:= Minimum(a,d)^(k-1)*PhiAL(N,a,d) / N^(k div 2-1);
    end if;
  end for;
  ret := -S1 / 2 - S2 / 2;
  if k eq 2 then
     ret +:= &+[n div d : d in Divisors(n) | GCD(d,N) eq 1];
  end if;
  return ret;
end function;

function TraceFormulaGamma0ALTrivialNew(N, k)
    ms := [d : d in Divisors(N) | N mod d^2 eq 0];
    trace := &+[Integers() | MoebiusMu(m)*TraceFormulaGamma0AL(1, N div m^2, k) : m in ms];
    return trace;
end function;

// At the moment only works for Hecke operators at primes
function TrivialContribution(N, k, p)
    assert IsPrime(p);
    N_primes_2 := [d : d in Divisors(N) | IsSquare(N*p div d) and ((N div d) mod p^3 ne 0) and (d mod p ne 0)];
    // trace_2 := &+[Integers() | get_trace(N_prime, k, 1 : New) : N_prime in N_primes_2];
    trace_2 := &+[Integers() | TraceFormulaGamma0ALTrivialNew(N_prime, k) : N_prime in N_primes_2];
    trace_3 := 0;
    if (N mod p eq 0) then
	N_primes_3 := [d : d in Divisors(N div p) | IsSquare(N div (d*p))];
	trace_3 := &+[Integers() | TraceFormulaGamma0ALTrivialNew(N_prime, k) :  N_prime in N_primes_3];
    end if;
    return p^(k div 2) * trace_3  - p^(k div 2 - 1)*trace_2;
end function;

// At the moment only works for Hecke operators at primes
function TraceFormulaGamma0ALNew(p, N, k)
    if (p eq 1) then return TraceFormulaGamma0ALTrivialNew(N, k); end if;
    assert IsPrime(p);
    ms := [d : d in Divisors(N) | (N mod d^2 eq 0) and (d mod p ne 0)];
    trace := &+[Integers() | MoebiusMu(m)*(TraceFormulaGamma0AL(p, N div m^2, k) - TrivialContribution(N div m^2, k, p)) : m in ms];
    return trace;
end function;

function A1(n,N,k)
  if (not IsSquare(n)) or (GCD(n,N) ne 1) then
    return 0;
  end if;
  return n^(k div 2 - 1)*phi1(N)*(k-1)/12;
end function;

function phi1(N)
    return N * &*[ Rationals() | 1 + 1/p : p in PrimeDivisors(N)];
end function;

function mu(N,t,f,n)
  N_f := GCD(N,f);
  primes := [x[1] : x in Factorization(N) | (N div N_f) mod x[1] ne 0];
  s := #[x : x in [0..N-1] | (GCD(x,N) eq 1) and ((x^2 - t*x+n) mod (GCD(f*N, N^2)) eq 0)];
  prod := IsEmpty(primes) select 1 else &*[ 1 + 1/p : p in  primes];
  assert N_f * prod eq (phi1(N) / phi1(N div GCD(N,f)));
  return N_f * prod * s;
end function;

function A2(n,N,k)
  R<x> := PolynomialRing(Rationals());
  max_abst := Floor(SquareRoot(4*n));
  if IsSquare(n) then max_abst -:= 1; end if;
  ret := 0;
  for t in [-max_abst..max_abst] do
      // print "t = ", t;
      F<rho> := NumberField(x^2 - t*x+n);
      rho_bar := t - rho;
      p := Rationals()!((rho^(k-1) - rho_bar^(k-1)) / (rho - rho_bar));
      assert p eq P(k,t,n);
      t_sum := 0;
      for d in Divisors(4*n - t^2) do
	  is_sq, f := IsSquare(d);
	  if is_sq then
              D := (t^2-4*n) div f^2;
	      if D mod 4 in [0,1] then
		  O := QuadraticOrder(BinaryQuadraticForms(D));
		  h := #PicardGroup(O);
		  w := #TorsionSubgroup(UnitGroup(O));
		  t_sum +:= (h/w) * mu(N,t,f,n);
	      end if;
	  end if;
      end for;
      // print "t_sum = ", t_sum;
      // print "p = ", p;
      ret -:= p*t_sum;
      // print "ret = ", ret;
  end for;
  return ret;
end function;

function A3(n,N,k)
  g := 1;
  ds := [d : d in Divisors(n) | d^2 lt n];
  ret := 0;
  for d in ds do
      // print "d = ", d;
      cs := [c : c in Divisors(N) |
	     (GCD(N div g, n div d - d) mod GCD(c, N div c) eq 0)
	    and (GCD(d mod c, c) eq 1) and (GCD((n div d) mod (N div c), (N div c)) eq 1)];
      ret -:= d^(k-1) * &+[Integers() | EulerPhi(GCD(c, N div c)) : c in cs];
      // print "ret = ", ret;
  end for;
  is_sq, d := IsSquare(n);
  if is_sq then
      cs := [c : c in Divisors(N) |
	     (GCD(N div g, n div d - d) mod GCD(c, N div c) eq 0)
	     and (GCD(d mod c, c) eq 1) and (GCD((n div d) mod (N div c), (N div c)) eq 1)];
      ret -:= 1/2 * d^(k-1) * &+[Integers() | EulerPhi(GCD(c, N div c)) : c in cs];
  end if;
  return ret;
end function;

function A4(n,N,k)
  if k eq 2 then
    return &+[t : t in Divisors(n) | GCD(N, n div t) eq 1];
  end if;
  return 0;
end function;

function TraceCohen(n,N,k)
  return A1(n,N,k) + A2(n,N,k) + A3(n,N,k) + A4(n,N,k);
end function;

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

// Formula from Popa
function TraceFormulaGamma0HeckeAL(N, k, n, Q)
    assert k ge 2;
    if (n eq 0) then return 0; end if; // for compatibility with q-expansions
    S1 := 0;
    Q_prime := N div Q;
    assert GCD(Q, Q_prime) eq 1;
    w := k - 2;
    max_abst := Floor(SquareRoot(4*Q*n)) div Q;
    for tQ in [-max_abst..max_abst] do
	t := tQ*Q;
	for u in Divisors(Q) do
	    for u_prime in Divisors(Q_prime) do
		if ((4*n*Q-t^2) mod (u*u_prime)^2 eq 0) then
		    // print "u =", u, " u_prime = ", u_prime, "t = ", t;
		    S1 +:= P(k,t,Q*n)*H((4*Q*n-t^2) div (u*u_prime)^2)*Cfast(Q_prime,u_prime,t,Q*n)
			   *MoebiusMu(u) / Q^(w div 2);
		    // print "S1 = ", S1;
		end if;
	    end for;
	end for;
    end for;
    S2 := 0;
    for d in Divisors(n*Q) do
	a := n*Q div d;
	if (a+d) mod Q eq 0 then
	    // print "a = ", a, "d = ", d;
	    S2 +:= Minimum(a,d)^(k-1)*Phil(N,Q,a,d) / Q^(w div 2);
	    // print "S2 = ", S2;
	end if;
    end for;
    // print "S2 = ", S2;
    ret := -S1 / 2 - S2 / 2;
    if k eq 2 then
	ret +:= &+[n div d : d in Divisors(n) | GCD(d,N) eq 1];
    end if;
    return ret;
end function;

function get_trace_hecke_AL(N, k, n, Q : New := false)
    C := CuspidalSubspace(ModularSymbols(N,k,1));
    if New then
	C := NewSubspace(C);
    end if;
    al := AtkinLehner(C,Q);
    T := HeckeOperator(C,n);
    return Trace(T*al);
end function;

function d_prime(d, N, Q, N_prime)
    return GCD(d, N div Q) * GCD(N div (N_prime*d), Q);
end function;

function Q_prime(N, Q, N_prime)
    return GCD(Q, N_prime);
end function;

function dd_prime(d, N, Q, N_prime, n)
    d_p := d_prime(d, N, Q, N_prime);
    if (GCD(d_p, n) * d) mod d_p ne 0 then
	return 0;
    end if;
    return (GCD(d_p, n) * d) div d_p;
end function;

function n_prime(d, N, Q, N_prime, n)
    d_p := d_prime(d, N, Q, N_prime);
    dd_p := dd_prime(d, N, Q, N_prime, n) * GCD(d_p, n);
    return n div dd_p;
end function;

function get_ds(N, Q, N_prime, n)
    divs := Divisors(N div N_prime);
    ret := [];
    for d in divs do
	d_p := d_prime(d, N, Q, N_prime);
	dd_p := dd_prime(d, N, Q, N_prime, n);
	if (dd_p eq 0) then
	    continue;
	end if;
	if (GCD(dd_p, N_prime) eq 1) and (n mod (dd_p * GCD(d_p, n)) eq 0) then
	    Append(~ret, d);
	end if;
    end for;
    return ret;
end function;

procedure testRelationNewSubspacesGeneral(N, k, n, Q)
    s := 0;
    for N_prime in Divisors(N) do
	ds := get_ds(N, Q, N_prime, n);
	for d in ds do
	    n_p := n_prime(d, N, Q, N_prime, n);
	    d_p := d_prime(d, N, Q, N_prime);
	    dd_p := dd_prime(d, N, Q, N_prime, n);
	    Q_p := Q_prime(N, Q, N_prime);
	    term := (n div n_p)^(k div 2 - 1);
	    term *:= GCD(d_p, n);
	    term *:= MoebiusMu(dd_p);
	    term *:= get_trace_hecke_AL(N_prime, k, n_p, Q_p : New);
	    s +:= term;
	end for;
    end for;
    assert s eq get_trace_hecke_AL(N, k, n, Q);
end procedure;
