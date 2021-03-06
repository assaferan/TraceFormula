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
// compare timings
// is_sq := KroneckerSymbol(y_0,p);
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
  return [x : x in [0..N-1] | GCD(x,N) eq 1 and (x^2 - t*x + n) mod N*u eq 0];
end function;

function phi1(N)
  primes := [f[1] : f in Factorization(N)];
  if IsEmpty(primes) then return N; end if;
  return Integers()!(N * &*[1 + 1/p : p in primes]);
end function;

function B(N, u, t, n)
// assert Sfast(N,u,t,n) eq #S(N,u,t,n);
//return #S(N,u,t,n) * phi1(N) div phi1(N div u);
  return Sfast(N,u,t,n) * phi1(N) div phi1(N div u);  
end function;

function C(N, M)
  a, b, c, d := Explode(M);
  G_M := GCD(c, d-a, b);
  return B(N, GCD(G_M, N), Trace(M), Determinant(M));
end function;

function C(N, u, t, n)
  return &+[B(N, u div d, t, n) * MoebiusMu(d) : d in Divisors(u)];
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

//ret := &+[ClassNumber(BinaryQuadraticForms(-n div d)) : d in Divisors(n) | IsSquare(d) and (n div d) mod 4 in [0,3] ];
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
      if (alpha ne 0) or (N eq 1) then
        ret +:= EulerPhi(g);
      end if;
    end if;
  end for;
  return ret;
end function;

// Gegenbauer polynomials
function P(k, t, m)
  R<x> := PowerSeriesRing(Rationals());
  return Coefficient((1 - t*x+m*x^2)^(-1), k-2);
end function;

// At the moment, this seems to work when N is prime
function TraceFormulaGamma0(n, N, k)
  S1 := 0;
  max_abst := Floor(SquareRoot(4*n));
  for t in [-max_abst..max_abst] do
    for u in Divisors(N) do
      if ((4*n-t^2) mod u^2 eq 0) then
	S1 +:= P(k,t,n)*H((4*n-t^2) div u^2)*C(N,u,t,n);
      end if;
    end for;
  end for;
  S2 := 0;
  for d in Divisors(n) do
    a := n div d;
    S2 +:= Minimum(a,d)^(k-1)*Phi(N,a,d);
  end for;
  ret := -S1 / 2 - S2 / 2;
  if k eq 2 then
     ret +:= &+[n div d : d in Divisors(n) | GCD(d,N) eq 1];
  end if;
  return ret;
end function;

function PhiAL(N, a, d)
  return EulerPhi(N) / N;
end function;

// This seems to work when N is prime, unless n eq N
function TraceFormulaGamma0AL(n, N, k)
  S1 := 0;
//  max_abst := Floor(SquareRoot(4*N*n));
  max_abst := Floor(SquareRoot(4*N*n)) div N;
  // ts := [t : t in [-max_abst..max_abst] | t mod N eq 0];
  // for t in ts do
  for tN in [-max_abst..max_abst] do
    t := tN*N;
    for u in Divisors(N) do
      if ((4*n-t^2) mod u^2 eq 0) then
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

function A1(n,N,k)
  if not IsSquare(n) then
    return 0;
  end if;
  return n^(k div 2 - 1)*phi1(N)*(k-1)/12;
end function;

function mu(N,t,f,n)
  N_f := GCD(N,f);
  primes := [x[1] : x in Factorization(N) | (N div N_f) mod x[1] ne 0];
  s := #[x : x in [0..N-1] | (x^2 - t*x+n) mod (GCD(f*N, N^2)) eq 0];
  prod := IsEmpty(primes) select 1 else &*[ 1 + 1/p : p in  primes];
  return N_f * prod * s;
end function;

// This is not working yet. Probably due to class number issues.
function A2(n,N,k)
  R<x> := PolynomialRing(Rationals());
  max_abst := Floor(SquareRoot(4*n));
  if IsSquare(n) then max_abst -:= 1; end if;
  ret := 0;
  for t in [-max_abst..max_abst] do
    F<rho> := NumberField(x^2 - t*x+n);
    rho_bar := t - rho;
    p := Rationals()!((rho^(k-1) - rho_bar^(k-1)) / (rho - rho_bar));
    assert p eq P(k,t,n);
    t_sum := 0;
    for d in Divisors(4*n - t^2) do
      is_sq, f := IsSquare(d);
      if is_sq then
        D := (t^2-4*n) div f^2;
        h := (D mod 4 in [0,1]) select ClassNumber(D) else 0;
        w := #UnitGroup(Integers(QuadraticField(D)));
        t_sum +:= (h/w) * mu(N,t,f,n);
      end if;
    end for;
    ret -:= p*t_sum;
  end for;
  return ret;
end function;

function A3(n,N,k)
  g := 1;
  ds := [d : d in Divisors(N) | d^2 lt n];
  ret := 0;
  for d in ds do  
    cs := [c : c in Divisors(N) |
	     GCD(N div g, n div d - d) mod GCD(c, N div c) eq 0];
    ret -:= d^(k-1) * &+[EulerPhi(GCD(c, N div c)) : c in cs];
  end for;
  is_sq, d := IsSquare(n);
  if is_sq then
    cs := [c : c in Divisors(N) |
	     GCD(N div g, n div d - d) mod GCD(c, N div c) eq 0];
    ret -:= 1/2 * d^(k-1) * &+[EulerPhi(GCD(c, N div c)) : c in cs];
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
