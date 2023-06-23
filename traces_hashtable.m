function Sfast_prime1(p,t,n, y)
  if p eq 2 then
    if IsOdd(1-t+n) then
      return 0;
    else
      return 1;
    end if;
  end if;
//  y := t^2-4*n;
  return KroneckerSymbol(y,p) + 1;
end function;

function Sfast_prime(p,u,t,n)
  if p eq 2 then
    if IsOdd(1-t+n) then
      return 0;
    elif u eq 1 then
      return 1;
    elif (1-t+n) mod 4 eq 0 then
      return IsEven(t) select 2 else 1;
    else
      return 0;
    end if;
  else
    y := t^2-4*n;
    if (u eq 1) or (y mod p ne 0) then
      return KroneckerSymbol(y,p) + 1;
    end if;
    if (y mod p^2) eq 0 then
      return 1;
    else
      return 0;
    end if;
  end if;
end function;

function Sfast(N, u, t, n)
// fac := Factorization(N*u);
  primes := PrimeDivisors(N);
  num_sols := 1;
//for f in fac do
  for p in primes do
    // p,e := Explode(f);
    if p eq 2 then
      // we try to postpone calculation of valuation
      if IsOdd(1-t+n) then
        return 0;
      end if;
      e2 := Valuation(1-t+n, 2);
      e := Valuation(N,p) + Valuation(u,p);
      if (e2 eq 0) or (e gt e2) then
        return 0;
      end if;
      if IsEven(t) then
        num_sols *:= 2^(e-1);
      end if;
    else
      y := t^2-4*n;
      // we try to postpone calculation of valuation
      // so as e ge 1 we handle e_y = 0
      if y mod p ne 0 then
        is_sq := KroneckerSymbol(y,p);
        if is_sq eq -1 then
	  return 0;
        end if;
        num_sols *:= 2;
      else
        e_y := Valuation(y, p);
        e := Valuation(N,p) + Valuation(u,p);
        if e_y lt e then
          if IsOdd(e_y) then
            return 0;
          end if;
	  y_0 := y div p^e_y;
          // This is faster
          is_sq := KroneckerSymbol(y_0,p);
          //is_sq := IsSquare(Integers(p)!y_0);
          //if not is_sq then
          if is_sq eq -1 then
	    return 0;
          end if;
          num_sols *:= p^(e_y div 2) * 2;
        else
	  num_sols *:= p^(e div 2);  
        end if;
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
 primes := PrimeDivisors(N);
 if IsEmpty(primes) then return N; end if;
 // This seems to be slightly faster
 return Integers()!(N * &*[1 + 1/p : p in primes]);
// prod := &* primes;
// mult := &*[1+p : p in primes];
// return (N div prod) * mult;
end function;

function phi1_quot(N, u)
  primes := [p : p in PrimeDivisors(N) | (N mod (p*u)) ne 0] ;
  if IsEmpty(primes) then return u; end if;
  return Integers()!(u* &*[1 + 1/p : p in primes]);
end function;

function B(N, u, t, n)
// assert Sfast(N,u,t,n) eq #S(N,u,t,n);
//return #S(N,u,t,n) * phi1(N) div phi1(N div u);
//  return Sfast(N,u,t,n) * phi1(N) div phi1(N div u);
//  tmp := phi1(N) div phi1(N div u);
//  tt := phi1_quot(N,u);
//  assert tmp eq tt;
  return Sfast(N,u,t,n) * phi1_quot(N,u);
end function;

function C(N, M)
  a, b, c, d := Explode(M);
  G_M := GCD(c, d-a, b);
  return B(N, GCD(G_M, N), Trace(M), Determinant(M));
end function;

function C(N, u, t, n)
  return &+[B(N, u div d, t, n) * MoebiusMu(d) : d in Divisors(u)];
end function;

function CC(N, u, D)
  fac := Factorization(N);
  res := 1;
  for f in fac do
    p, a := Explode(f);
    i := Valuation(u, p);
    if i eq 0 then continue; end if;
    if i eq a then
      res *:= p^((a+1) div 2);
      continue;
    end if;
    b := Valuation(D,p);
    if p eq 2 then
      if (i le b-a-2) and (i ge 1) and IsEven(i-a) then
	res *:= 2^((i-1) div 2);
      elif (i eq b-a-1) and IsEven(i-a) then
        res *:= -2^((i-1) div 2);
      elif (i eq b-a) and IsEven(i-a) then
        res *:= 2^((i-1) div 2) * KroneckerCharacter(-4)(D div 2^b);
      elif (i eq b-a+1) and IsOdd(i-a) and ((D div 2^b) mod 4 eq 1) then
        res *:= 2^(i div 2) * KroneckerSymbol(D div 2^b, 2);
      else
        return 0;
      end if;
    else
      if (i le b-a) and (i ge 1) and IsEven(i-a) then
	res *:= p^((i+1) div 2) - p^((i-1) div 2);
      elif (i eq b-a+1) and IsEven(i-a) then
        res *:= -p^((i-1) div 2);
      elif (i eq b-a+1) and IsOdd(i-a) then
        res *:= p^(i div 2) *KroneckerSymbol(D div p^b, p);
      else
	return 0;
      end if;
    end if;
  end for;
  return res;
end function;

// version for primes
function CC_prime(p, u, D)
  return u;
end function;

function Cfaster(N, u, t, n, s)
  if s eq 0 then
    return 0;
  end if;
//  assert CC_prime(N, u, t^2-4*n) eq CC(N, u, t^2-4*n);
//  return s * CC(N, u, t^2-4*n);
  return s*u;
end function;

function Cfast(N, u, t, n)
// assert Sfast(N,1,t,n) eq Sfast_prime1(N,t,n);
// s_fast := Sfast(N,1,t,n);
  s_fast := Sfast_prime1(N,t,n);
  if s_fast eq 0 then
     return 0;
  end if;
  return s_fast * CC(N, u, t^2-4*n);
end function;

procedure H(n, ~stored_values, ~ret)
  is_in := IsDefined(stored_values, n);
  if is_in then
    ret := stored_values[n];
  else
    if n lt 0 then
      is_sq, u := IsSquare(-n);
      ret := (is_sq select -u/2 else 0);
      stored_values[n] := ret;
    end if;
    if n eq 0 then
      ret := -1/12;
      stored_values[n] := ret;
    end if;
    if (not IsDefined(stored_values, n)) and (n mod 4 in [1,2]) then
      ret := 0;
      stored_values[n] := ret;
    end if;

    if not IsDefined(stored_values, n) then
    //ret := &+[ClassNumber(BinaryQuadraticForms(-n div d)) : d in Divisors(n) | IsSquare(d) and (n div d) mod 4 in [0,3] ];
      ret := &+[ClassNumber(-n div d) : d in Divisors(n)
    		       | IsSquare(d) and (n div d) mod 4 in [0,3] ];
      if IsSquare(n) and IsEven(n) then
        ret -:= 1/2;
      end if;
      if n mod 3 eq 0 and IsSquare(n div 3) then
        ret -:= 2/3;
      end if;
      stored_values[n] := ret;
    end if;
  end if;
end procedure;

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
  if (k eq 2) then
    return 1;
  end if;
  R<x> := PowerSeriesRing(Rationals());
  return Coefficient((1 - t*x+m*x^2)^(-1), k-2);
end function;

// At the moment, this seems to work when N is prime
procedure TraceFormulaGamma0(n, N, k, ~stored_H_values, ~ret)
  Hval := false;
  S1 := 0;
  max_abst := Floor(SquareRoot(4*n));
//  for t in [-max_abst..max_abst] do
  for t in [0..max_abst] do
    scalar := (t eq 0) select 1 else 2;
    t_sum := 0;
    y := 4*n-t^2;
    s := Sfast_prime1(N,t,n,-y);
    if s ne 0 then
      for u in Divisors(N) do
        if (y mod u^2 eq 0) then
          H(y div u^2, ~stored_H_values, ~Hval);
          cfast := Cfaster(N,u,t,n,s);
// S1 +:= scalar*P(k,t,n)*Hval*cfast;
          t_sum +:= Hval*cfast;
        end if;
      end for;
      S1 +:= scalar*t_sum;
    end if;
  end for;
  S2 := 0;
  for d in Divisors(n) do
    a := n div d;
// phi := Phi(N,a,d);
//    assert (phi eq 2) or GCD(n,N) gt 1;
//    S2 +:= Minimum(a,d)^(k-1)*Phi(N,a,d);
    S2 +:= Minimum(a,d)*2;
  end for;
  ret := -S1 / 2 - S2 / 2;
  if k eq 2 then
     ret +:= &+[n div d : d in Divisors(n) | GCD(d,N) eq 1];
  end if;
  //return ret;
end procedure;

function PhiAL(N, a, d)
  return EulerPhi(N) / N;
end function;

// This seems to work when N is prime, unless n eq N
procedure TraceFormulaGamma0AL(n, N, k, ~stored_H_values, ~ret)
  Hval := false;
  S1 := 0;
//  max_abst := Floor(SquareRoot(4*N*n));
  max_abst := Floor(SquareRoot(4*N*n)) div N;
  // ts := [t : t in [-max_abst..max_abst] | t mod N eq 0];
  // for t in ts do
  for tN in [-max_abst..max_abst] do
    t := tN*N;
    for u in Divisors(N) do
      if ((4*N*n-t^2) mod u^2 eq 0) then
        H((4*N*n-t^2) div u^2, ~stored_H_values, ~Hval);
        S1 +:= P(k,t,N*n)* Hval *C(1,1,t,N*n) * MoebiusMu(u) / N^(k div 2-1);
      end if;
    end for;
  end for;

// print "S1 = ", S1;
  S2 := 0;
  for d in Divisors(n*N) do
    a := n*N div d;
    if (a+d) mod N eq 0 then 
      S2 +:= Minimum(a,d)^(k-1)*PhiAL(N,a,d) / N^(k div 2-1);
    end if;
  end for;
// print "S2 = ", S2;
  ret := -S1 / 2 - S2 / 2;
  if k eq 2 then
     ret +:= &+[n div d : d in Divisors(n) | GCD(d,N) eq 1];
  end if;
  //return ret;
end procedure;

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

procedure ComputeTraces(p, ~stored_values, ~res)
  for n in [1..#res] do
    TraceFormulaGamma0(n,p,2,~stored_values,~res[n]);
  end for;
end procedure;

procedure ComputeTracesAL(p, ~stored_values, ~res)
  for n in [1..#res] do
    TraceFormulaGamma0AL(n,p,2,~stored_values,~res[n]);
  end for;
end procedure;

function GetTraces(p)
  stored_values := AssociativeArray();
  res := [0 : n in [1..p div 12]];
  ComputeTraces(p, ~stored_values, ~res);
  return res;
end function;

function GetTracesAL(p)
  stored_values := AssociativeArray();
  res := [0 : n in [1..p div 12]];
  ComputeTracesAL(p, ~stored_values, ~res);
  return res;
end function;

procedure testTraces(indices, values)
  stored_values := AssociativeArray();
  traces_full := [Rationals() | 0 : n in [1..1000]];
  traces_AL := [Rationals() | 0 : n in [1..1000]];
  for j in [1..#indices] do
    idx := indices[j];
    p := idx[1];
    sign := idx[2];
    ComputeTraces(p, ~stored_values, ~traces_full);
    ComputeTracesAL(p, ~stored_values, ~traces_AL);
    traces := [(traces_full[i] + sign*traces_AL[i]) / 2 : i in [1..1000]];
    assert &and[traces[i] eq values[j][i] : i in [1..1000] | i mod p ne 0];
    if j mod 100 eq 0 then
      print j, "/", #indices;
    end if;
  end for;
end procedure;

procedure testPariTraces(indices, values)
  for j in [1..#indices] do
    idx := indices[j];
    p := idx[1];
    sign := idx[2];
    traces_full := eval(Sprintf("traces_%o", p));
    traces_AL := eval(Sprintf("tracesAL_%o", p));
    traces := [(traces_full[i] + sign*traces_AL[i]) / 2 : i in [2..1000]];
    assert &and[traces[i] eq values[j][i] : i in [1..1000] | i mod p ne 0];
    if j mod 100 eq 0 then
      print j, "/", #indices;
    end if;
  end for;
end procedure;

function timeTraces(upTo)
  tt := Cputime();
  stored_values := AssociativeArray();
  for p in PrimesUpTo(upTo) do
    traces_full := [0 : n in [1..p div 12]];
    traces_AL := [0 : n in [1..p div 12]];
    ComputeTraces(p, ~stored_values, ~traces_full);
    ComputeTracesAL(p, ~stored_values, ~traces_AL);
  end for;
  return Cputime() - tt;
end function;
