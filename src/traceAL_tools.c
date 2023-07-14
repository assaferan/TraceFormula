#include <pari/pari.h>
#include <inttypes.h>
#include <time.h>

#include "traceAL_tools.h"


/*
GP;install("traceAL", "GG&&", "traceAL", "./libtraceal.so");
*/

/* fa = factorization of -D > 0, return -D0 > 0 (where D0 is fundamental) */
static long corediscs_fact(GEN fa)
{
  GEN P = gel(fa,1), E = gel(fa,2);
  long i, l = lg(P), m = 1;
  for (i = 1; i < l; i++)
  {
    long p = P[i], e = E[i];
    if (e & 1) m *= p;
  }
  if ((m&3L) != 3) m <<= 2;
  return m;
}

static void constfact(long lim);
static void constdiv(long lim);
static void consttabh(long lim);
static void constcoredisc(long lim);

static THREAD cache caches[4] = {
{ "Factors",  NULL,  50000,    50000, &constfact, 0, 0, 0 },
{ "Divisors", NULL,  50000,    50000, &constdiv, 0, 0, 0 },
{ "H",        NULL, 100000, 10000000, &consttabh, 0, 0, 1 },
{ "CorediscF",NULL, 100000, 10000000, &constcoredisc, 0, 0, 0 }
};

static void cache_reset(long id) { caches[id].miss = caches[id].maxmiss = 0; }
static void cache_delete(long id) {/*if (caches[id].cache != NULL) gunclone(caches[id].cache); */ }
static void cache_set(long id, GEN S)
{
  GEN old = caches[id].cache;
  caches[id].cache = gclone(S);
  //  guncloneNULL(old);
  if (old != NULL)
    gunclone(old);
}

/* handle a cache miss: store stats, possibly reset table; return value
 * if (now) cached; return NULL on failure. HACK: some caches contain an
 * W64 where the 0 value is impossible, and return it (typecast to GEN) */
static GEN cache_get(long id, W64 D)
{
#ifdef DEBUG
  printf("In cache_get, with id = %ld, D = %" PRId64 "\n", id, D);
#endif // DEBUG
  cache *S = &caches[id];
  const W64 d = S->compressed? D>>1: D;
  W64 max, l;

   if (!S->cache)
  {
    max = maxuu(minuu(D, S->maxself), S->minself);
    S->init(max);
    l = lg(S->cache);
  }
  else
  {
    l = lg(S->cache);
    if (l <= d)
    {
      if (D > S->maxmiss) S->maxmiss = D;
      if (DEBUGLEVEL >= 3)
        err_printf("miss in cache %s: %lu, max = %lu\n",
                   S->name, D, S->maxmiss);
      if (S->miss++ >= 5 && D < S->maxself)
      {
        max = minuu(S->maxself, (long)(S->maxmiss * 1.2));
        if (max <= S->maxself)
        {
          if (DEBUGLEVEL >= 3)
            err_printf("resetting cache %s to %lu\n", S->name, max);
	  S->init(max); l = lg(S->cache);
        }
      }
    }
  }
  return (l <= d)? NULL: gel(S->cache, d);
}
static GEN cache_report(long id)
{
  cache *S = &caches[id];
  GEN v = zerocol(5);
  gel(v,1) = strtoGENstr(S->name);
  if (S->cache)
  {
    gel(v,2) = utoi(lg(S->cache)-1);
    gel(v,3) = utoi(S->miss);
    gel(v,4) = utoi(S->maxmiss);
    gel(v,5) = utoi(gsizebyte(S->cache));
  }
  return v;
}

GEN getcache(void)
  {
  pari_sp av = avma;
  GEN M = cgetg(6, t_MAT);
  gel(M,1) = cache_report(cache_FACT);
  gel(M,2) = cache_report(cache_DIV);
  gel(M,3) = cache_report(cache_H);
  gel(M,4) = cache_report(cache_D);
  return gerepilecopy(av, shallowtrans(M));
}

void pari_close_mf(void)
{
  cache_delete(cache_FACT);
  cache_delete(cache_DIV);
  cache_delete(cache_H);
  cache_delete(cache_D);
}

/*************************************************************************/
/* a odd, update local cache (recycle memory) */
static GEN update_factor_cache(long a, long lim, long *pb)
{
  const long step = 16000; /* even; don't increase this: RAM cache thrashing */
  if (a + 2*step > lim)
    *pb = lim; /* fuse last 2 chunks */
  else
    *pb = a + step;
  return vecfactoroddu_i(a, *pb);
}

/* assume lim < MAX_LONG/8 */
static void constcoredisc(long lim)
{
  pari_sp av2, av = avma;
  GEN D = caches[cache_D].cache, CACHE = NULL;
  long cachea, cacheb, N, LIM = !D ? 4 : lg(D)-1;
  if (lim <= 0) lim = 5;
  if (lim <= LIM) return;
  cache_reset(cache_D);
  D = zero_zv(lim);
  av2 = avma;
  cachea = cacheb = 0;
  for (N = 1; N <= lim; N+=2)
  { /* N odd */
    long i, d, d2;
    GEN F;
    if (N > cacheb)
      {
	//set_avma(av2);
	avma = av2;
	cachea = N;
      CACHE = update_factor_cache(N, lim, &cacheb);
    }
    F = gel(CACHE, ((N-cachea)>>1)+1); /* factoru(N) */
    D[N] = d = corediscs_fact(F); /* = 3 mod 4 or 4 mod 16 */
    d2 = odd(d)? d<<3: d<<1;
    for (i = 1;;)
    {
      if ((N << i) > lim) break;
      D[N<<i] = d2; i++;
      if ((N << i) > lim) break;
      D[N<<i] = d; i++;
    }
  }
  cache_set(cache_D, D);
  avma = av; // set_avma(av);
}

static void constfact(long lim)
{
  pari_sp av;
  GEN VFACT = caches[cache_FACT].cache;
  long LIM = VFACT? lg(VFACT)-1: 4;
  if (lim <= 0) lim = 5;
  if (lim <= LIM) return;
  cache_reset(cache_FACT); av = avma;
  cache_set(cache_FACT, vecfactoru_i(1,lim)); avma = av; // set_avma(av);
}

static void constdiv(long lim)
{
  pari_sp av;
  GEN VFACT, VDIV = caches[cache_DIV].cache;
  long N, LIM = VDIV? lg(VDIV)-1: 4;
  if (lim <= 0) lim = 5;
  if (lim <= LIM) return;
  constfact(lim);
  VFACT = caches[cache_FACT].cache;
  cache_reset(cache_DIV); av = avma;
  VDIV  = cgetg(lim+1, t_VEC);
  for (N = 1; N <= lim; N++) gel(VDIV,N) = divisorsu_fact(gel(VFACT,N));
  cache_set(cache_DIV, VDIV); avma = av; // set_avma(av);
}

/* n > 1, D = divisors(n); sets L = 2*lambda(n), S = sigma(n) */
static void lamsig(GEN D, long *pL, long *pS)
{
  pari_sp av = avma;
  long i, l = lg(D), L = 1, S = D[l-1]+1;
  for (i = 2; i < l; i++) /* skip d = 1 */
  {
    long d = D[i], nd = D[l-i]; /* nd = n/d */
    if (d < nd) { L += d; S += d + nd; }
    else
    {
      L <<= 1; if (d == nd) { L += d; S += d; }
      break;
    }
  }
  //set_avma(av);
  avma = av;
  *pL = L; *pS = S;
}

/* table of 6 * Hurwitz class numbers D <= lim */
static void consttabh(long lim)
{
  pari_sp av = avma, av2;
  GEN VHDH0, VDIV, CACHE = NULL;
  GEN VHDH = caches[cache_H].cache;
  long r, N, cachea, cacheb, lim0 = VHDH? lg(VHDH)-1: 2, LIM = lim0 << 1;

  if (lim <= 0) lim = 5;
  if (lim <= LIM) return;
  cache_reset(cache_H);
  r = lim&3L; if (r) lim += 4-r;
  cache_get(cache_DIV, lim);
  VDIV = caches[cache_DIV].cache;
  VHDH0 = cgetg(lim/2 + 1, t_VECSMALL);
  VHDH0[1] = 2;
  VHDH0[2] = 3;
  for (N = 3; N <= lim0; N++) VHDH0[N] = VHDH[N];
  av2 = avma;
  cachea = cacheb = 0;
    for (N = LIM + 3; N <= lim; N += 4)
  {
    long s = 0, limt = usqrt(N>>2), flsq = 0, ind, t, L, S;
    GEN DN, DN2;
    if (N + 2 >= lg(VDIV))
    { /* use local cache */
      GEN F;
      if (N + 2 > cacheb)
      {
        // set_avma(av2);
	avma = av2;
	cachea = N;
        CACHE = update_factor_cache(N, lim+2, &cacheb);
      }
      F = gel(CACHE, ((N-cachea)>>1)+1); /* factoru(N) */
      DN = divisorsu_fact(F);
      F = gel(CACHE, ((N-cachea)>>1)+2); /* factoru(N+2) */
      DN2 = divisorsu_fact(F);
    }
    else
    { /* use global cache */
      DN = gel(VDIV,N);
      DN2 = gel(VDIV,N+2);
    }
        ind = N >> 1;
    for (t = 1; t <= limt; t++)
    {
      ind -= (t<<2)-2; /* N/2 - 2t^2 */
      if (ind) s += VHDH0[ind]; else flsq = 1;
    }
    lamsig(DN, &L,&S);
    VHDH0[N >> 1] = 2*S - 3*L - 2*s + flsq;
    s = 0; flsq = 0; limt = (usqrt(N+2) - 1) >> 1;
    ind = (N+1) >> 1;
       for (t = 1; t <= limt; t++)
    {
      ind -= t<<2; /* (N+1)/2 - 2t(t+1) */
      if (ind) s += VHDH0[ind]; else flsq = 1;
    }
    lamsig(DN2, &L,&S);
    VHDH0[(N+1) >> 1] = S - 3*(L >> 1) - s - flsq;
  }
    cache_set(cache_H, VHDH0); avma = av; // set_avma(av);
}

/*************************************************************************/
/* Core functions using factorizations, divisors of class numbers caches */
static GEN myfactoru(long N)
{
  GEN z = cache_get(cache_FACT, N);
  return z? gcopy(z): factoru(N);
}

/* write n = mf^2. Return m, set f. */
static W64 mycore(W64 n, long *pf)
{
  pari_sp av = avma;
  GEN fa = myfactoru(n), P = gel(fa,1), E = gel(fa,2);
  long i, l = lg(P), m = 1, f = 1;
  for (i = 1; i < l; i++)
  {
    long j, p = P[i], e = E[i];
    if (e & 1) m *= p;
    for (j = 2; j <= e; j+=2) f *= p;
  }
  *pf = f; avma = av; return m; // return gc_long(av,m);
}

/* write -n = Df^2, D < 0 fundamental discriminant. Return D, set f. */
static Z64 mycoredisc2neg(W64 n, long *pf)
{
  W64 m, D = (W64)cache_get(cache_D, n);
  if (D) { *pf = usqrt(n/D); return -(Z64)D; }
  m = mycore(n, pf);
  if ((m&3) != 3) { m <<= 2; *pf >>= 1; }
  return (Z64)-m;
}

/* 1+p+...+p^e, e >= 1 */
static W64 usumpow(W64 p, Z64 e)
{
  W64 q = 1+p;
  Z64 i;
  for (i = 1; i < e; i++) q = p*q + 1;
  return q;
}

/* Hurwitz(D0 F^2)/ Hurwitz(D0)
 * = \sum_{f|F}  f \prod_{p|f} (1-kro(D0/p)/p)
 * = \prod_{p^e || F} (1 + (p^e-1) / (p-1) * (p-kro(D0/p))) */
static long get_sh(long F, long D0)
{
  GEN fa = myfactoru(F), P = gel(fa,1), E = gel(fa,2);
  long i, l = lg(P), t = 1;
  for (i = 1; i < l; i++)
  {
    long p = P[i], e = E[i], s = kross(D0,p);
    if (e == 1) { t *= 1 + p - s; continue; }
    if (s == 1) { t *= upowuu(p,e); continue; }
    t *= 1 + usumpow(p,e-1)*(p-s);
  }
  return t;
}

/* d > 0, d = 0,3 (mod 4). Return 6*hclassno(d); -d must be fundamental
 * Faster than quadclassunit up to 5*10^5 or so */
static ulong hclassno6u_count(W64 d)
{
  W64 a, b, b2, h = 0;
  Z64 f = 0;

  if (d > 500000)
    return 6 * itou(gel(quadclassunit0(utoineg(d), 0, NULL, 0), 1));

  /* this part would work with -d non fundamental */
  b = d&1; b2 = (1+d)>>2;
  if (!b)
  {
    for (a=1; a*a<b2; a++)
      if (b2%a == 0) h++;
    f = (a*a==b2); b=2; b2=(4+d)>>2;
     }
  while (b2*3 < d)
  {
    if (b2%b == 0) h++;
    for (a=b+1; a*a < b2; a++)
      if (b2%a == 0) h += 2;
    if (a*a == b2) h++;
    b += 2; b2 = (b*b+d)>>2;
  }
  if (b2*3 == d) return 6*h+2;
  if (f) return 6*h+3;
  return 6*h;
}

/* D > 0; 6 * hclassno(D), using D = D0*F^2 */
static long hclassno6u_2(W64 D, Z64 D0, Z64 F)
{
  long h;
#ifdef DEBUG
  printf("In hclassno6u_2, with D = %" PRId64 ", D0 = %" PRId64 ", F = %" PRId64 "\n", D, D0, F);
#endif // DEBUG
  if (F == 1) h = hclassno6u_count(D);
  else
  { /* second chance */
    h = (W64)cache_get(cache_H, -D0);
    if (!h) h = hclassno6u_count(-D0);
    h *= get_sh(F,D0);
  }
  return h;
}

/* D > 0; 6 * hclassno(D) (6*Hurwitz). Beware, cached value for D (=0,3 mod 4)
 * is stored at D>>1 */
ulong hclassno6w64(W64 D)
{
#ifdef DEBUG
  printf("In hclassno6w64, with D = %" PRId64 "\n", D);
#endif // DEBUG
  W64 z = (W64)cache_get(cache_H, D);
#ifdef DEBUG
  printf("From cache got z = H(D) = %" PRId64 "\n", z);
#endif // DEBUG
  Z64 D0;
  long F;
  if (z) return z;
  D0 = mycoredisc2neg(D, &F);
#ifdef DEBUG
  printf("Fundamental disc = D0 = %" PRId64 ", F = %ld\n", D0, F);
#endif // DEBUG
  return hclassno6u_2(D,D0,F);
}

W64 gtoW64(GEN D)
{
  GEN M = powgi(mkintn(1,2), mkintn(1,32));
  GEN D_L = gmod(D,M);
  GEN D_H = gdivent(D,M);
  
  W64 ret = gtou(D_H);
  ret <<= 32;
  ret |= gtou(D_L);
  return ret;
}

// returns 12*H
long H12(GEN D)
{
  // !! TODO - There might be an overflow here - try to fix that!
  GEN uu;
  long u;
  if (D == gen_0) return -1;
  if (D > gen_0)
    switch (gtos(gmod(D,mkintn(1,4)))) {
    case 0:
    case 3:
      return 2*hclassno6w64(gtoW64(D));
    case 1:
    case 2:
      return 0;
    }
  long is_sq = issquareall(gneg(D), &uu);
  u = gtos(uu);
  return (is_sq ? -6*u : 0); 
}

GEN mksintn(long l, long x)
{
  GEN res = gen_0;
   
  if (x > 0)
    res = gadd(res, mkintn(1,x));
  else
    res = gsub(res, mkintn(1,-x));

  return res;
}

GEN polyGegenbauer(long k, GEN t, GEN m)
{
  GEN pol_one = mkpoln(1,gen_1);
  GEN pol_quad = mkpoln(3,m,gneg(t),gen_1);
  GEN inv_pol = mkrfrac(pol_one, pol_quad);
  GEN inv_pol_ser = Ser0(inv_pol, -1, mkintn(1,(k-2) + 1), (k-2) + 1);
#ifdef DEBUG
  pari_printf("poly Gegenbauer = %Ps\n", inv_pol);
  pari_printf("poly Gegenbauer = %Ps\n", inv_pol_ser);
#endif // DEBUG
  GEN ret = truecoeff(inv_pol_ser, k-2);
#ifdef DEBUG
  pari_printf("poly Gegenbauer coeff = %Ps\n", ret);
#endif // DEBUG
  return ret;
}

/*
long alpha(ulong n)
{
  pari_sp av = avma;
  GEN fa = myfactoru(n), P = gel(fa,1), E = gel(fa,2);
  long i, l = lg(P), m = 1;
  for (i = 1; i < l; i++)
  {
    //     long j, e = E[i];
    long e = E[i];
    if ((e == 1) || (e == 2))
      m = -m;
    if (e >= 4)
      m = 0;
  }
  avma = av;
  return m;
}
*/

// In what follows we assume k >=2 is even

GEN traceAL(long N, long n, long k)
{
  GEN NN = mkintn(1,N);
  GEN nn = mkintn(1,n); 
  // const long nN = n*N;
  GEN nN = gmul(NN, nn); 
  //  const long n4N = nN << 2;
  GEN n4N = gmul(nN,mkintn(1,4));

  // Cohen does something wiser - see if it works here
  long limt, tN;
  GEN ret;
  GEN S1 = gen_0;

#ifdef DEBUG
  printf("In traceAL, k = %ld, N = %ld, n = %ld\n", k, N, n);
  pari_printf("n4N = %Ps, sqrt(n4N) = %Ps\n", n4N, gsqrt(n4N,3));
#endif // DEBUG
  limt = gtos(gdivent(gfloor(gsqrt(n4N,3)),NN));
#ifdef DEBUG
  printf("limt = %ld\n", limt);
#endif // DEBUG
  GEN div_nN = divisors(nN);
#ifdef DEBUG
  pari_printf("div_nN = %Ps\n", div_nN);
#endif // DEBUG
  GEN div_n = divisors(nn);
  GEN div_N = divisors(NN);
  long num_divs_nN = lg(div_nN);
  long num_divs_n = lg(div_n);
  long num_divs_N = lg(div_N);
  long phi = gtos(eulerphi(NN));
  GEN denom = powgi(NN, mkintn(1,(k/2)-1));
  for (tN = -limt ; tN <= limt; tN++) /* t^2 < 4Nn */
  {
    GEN t = gmul(NN, mksintn(1,tN));
#ifdef DEBUG
    pari_printf("tN = %ld, NN = %Ps, t = %Ps\n", tN, NN, t);
#endif // DEBUG
    GEN t2 = gmul(t,t);
    GEN D = gsub(n4N, t2);
    // pari_printf("t = %Ps, D = %Ps, ", t, D);
    GEN inner_sum_t = gen_0;
    for (long idx = 1; idx < num_divs_N; idx++) {
      GEN u = gel(div_N, idx);
      GEN u2 = gmul(u,u);
      if (gmod(D,u2) == gen_0) {
	inner_sum_t = gaddgs(inner_sum_t, moebius(u)*H12(gdivent(D,u2)));
      }
#ifdef DEBUG
      pari_printf("u = %Ps, H12(D / u^2) = %ld\n", u, H12(gdivent(D,u2)));
#endif // DEBUG
    }
    inner_sum_t = gmul(inner_sum_t, polyGegenbauer(k,t,nN));
    inner_sum_t = gdiv(inner_sum_t, denom);
    S1 = gadd(S1, inner_sum_t);
  }

#ifdef DEBUG
  pari_printf("Sum of class numbers is: %Ps\n", S1);
#endif // DEBUG
  GEN S2 = gen_0;

#ifdef DEBUG
  printf("num_divs_nN = %ld\n", num_divs_nN);
#endif // DEBUG
  for (long idx = 1; idx < num_divs_nN; idx++)
    {
#ifdef DEBUG
    pari_printf("div_nN[%ld] = %Ps\n", idx, gel(div_nN,idx));
#endif // DEBUG
    GEN d = gel(div_nN, idx);
    GEN a = gdivent(nN, d);
    if (gmod(gadd(a,d), NN) == gen_0)
    {
#ifdef DEBUG
      printf("Adding to S2...\n");
#endif // DEBUG
      S2 = gadd(S2, powgi(gmin(a,d), mkintn(1,k-1)));
    }
  }

#ifdef DEBUG
  pari_printf("Sum of divisors is: %Ps\n", S2);
#endif // DEBUG
  S2 = gmulgs(S2, 12*phi);
  S2 = gdivgs(S2, N);
  S2 = gdiv(S2, denom);

  ret = gadd(S1, S2);
  ret = gsub(gen_0, ret);
  ret = gdivgs(ret, 24);
  
  if (k == 2) {
    for (long idx = 1; idx < num_divs_n; idx++)
      {
	GEN d = gel(div_n, idx);
	if (ugcd(N,gtos(d)) == 1)
	  ret = gadd(ret, gdiv(nn,d));
      }
  }
  return ret;
}

GEN traceALNewTrivial(long N, long k)
{
  GEN trace = gen_0;
  GEN NN = mkintn(1,N);
  GEN div_N = divisors(NN);
  long num_divs_N = lg(div_N);

  for (long idx = 1; idx < num_divs_N; idx++) {
    GEN D = gel(div_N, idx);
    GEN D2 = gmul(D,D);
    
    // if ((N % (d*d)) != 0) continue;
    if (gmod(NN, D2) != gen_0) continue;
    long N_div_d2 = gtos(gdivent(NN, D2));
    trace = gadd(trace, gmulsg(moebius(D),traceAL(N_div_d2, 1, k)));
  }

#ifdef DEBUG
  pari_printf("Trace of W_{%Ps} on new subspace is: %Ps\n", NN, trace);
#endif // DEBUG
  
  return trace;
}

GEN traceALNewTrivialContribution(long N, long p, long k)
{
  GEN trace;
  GEN trace2 = gen_0;
  GEN trace3 = gen_0;
  GEN div_N = divisors(mkintn(1,N));
  long num_divs_N = lg(div_N);

  for (long idx = 1; idx < num_divs_N; idx++) {
    long d = gtos(gel(div_N, idx));
    long N_div_d = N / d;
    // We use valuation, since p^3 might not be a single-word integer
    if ((d % p == 0) || (gvaluation(mkintn(1,N_div_d), mkintn(1,p)) >= 3)) continue;
    if (! issquare(mkintn(1,p*N_div_d)) ) continue;
    trace2 = gadd(trace2, traceALNewTrivial(d, k));
  }
#ifdef DEBUG
  pari_printf("Contribution from N2 for level %d for prime %d and weight %d is: %Ps \n", N, p, k, trace2);
#endif // DEBUG

  if (N % p == 0) {
    long N_div_p = N / p;
    GEN div_N_div_p = divisors(mkintn(1,N_div_p));
    long num_divs_N_div_p = lg(div_N_div_p);
    for (long idx = 1; idx < num_divs_N_div_p; idx++) {
      long d = gtos(gel(div_N_div_p, idx));
      if (! issquare(mkintn(1,N_div_p / d)) ) continue;
      trace3 = gadd(trace3, traceALNewTrivial(d, k));
    }
  }
#ifdef DEBUG
  pari_printf("Contribution from N3 for level %d for prime %d and weight %d is: %Ps \n", N, p, k, trace3);
#endif // DEBUG
  
  // Since p is at most one word length, we can estimate the size of its powers accordingly
  trace = gsub(gmul(gpow(mkintn(1,p), mkintn(1,k/2), k/2), trace3),
	       gmul(gpow(mkintn(1,p), mkintn(1,k/2-1), k/2-1), trace2));

#ifdef DEBUG
  pari_printf("Trivial contribution from level %d for prime %d and weight %d is: %Ps \n", N, p, k, trace);
#endif // DEBUG
  return trace;
}

// At the moment only works for primes
GEN traceALNew(long N, long p, long k)
{
  GEN trace = gen_0;
  GEN NN = mkintn(1,N);
  GEN div_N = divisors(NN);
  long num_divs_N = lg(div_N);

#ifdef DEBUG
  printf("In traceALNew with N = %ld.\n", N);
#endif // DEBUG
  
  for (long idx = 1; idx < num_divs_N; idx++) {
    GEN D = gel(div_N, idx);
    GEN D2 = gmul(D,D);

#ifdef DEBUG
    pari_printf("D = %Ps \n", D);
    pari_printf("N mod (D^2) = %Ps \n", gmod(NN,D2));
    printf("(N mod (D^2) != 0) = %d \n", (gmod(NN,D2) != 0));
    printf("gen_0 (N mod (D^2) != 0) = %d \n", (gmod(NN,D2) != gen_0));
    printf("cmp (N mod (D^2) != 0) = %d \n", gcmp(gmod(NN,D2),gen_0));
    printf("gtos(D) mod p = %ld\n", gtos(D) % p);
    printf("(gtos(D) mod p == 0) =  %d\n", gtos(D) % p == 0);
#endif // DEBUG
    
    if (gmod(NN, D2) != gen_0) continue;
    if (gtos(D) % p == 0) continue;
    
    long N_div_d2 = gtos(gdivent(NN, D2));
    trace = gadd(trace, gmulsg(moebius(D), gsub(traceAL(N_div_d2,p,k), traceALNewTrivialContribution(N_div_d2,p,k))));
#ifdef DEBUG
    pari_printf("Accumulated trace is: %Ps \n", trace);
#endif // DEBUG
  }
  
  return trace;
}

GEN traceALprimes(long N, long k, long prec, int newspace, long start)
{
  GEN p_list = primes0(mkvec2(mkintn(1,start), nextprime(mkintn(1,prec))));
  long num_primes = lg(p_list);
  // adding also the trace for T_1 = 1
  // GEN res = cgetg(num_primes-1, t_VEC);
  GEN res = cgetg(num_primes, t_VEC);
  long p;

  GEN (*trace_func)(long, long, long);
  trace_func = (newspace ? &traceALNew : &traceAL);

#ifdef DEBUG
  printf("In traceALprimes. num_primes = %ld. newspace = %d. start = %ld. \n", num_primes, newspace, start);
#endif // DEBUG

  gel(res, 1) = (newspace ? traceALNewTrivial(N,k) : traceAL(N,1,k));
  
  for (long idx = 1; idx < num_primes - 1; idx++)
  {
    p = gtos(gel(p_list, idx));
    // gel(res, idx) = (*trace_func)(N, p, k);
    gel(res, idx+1) = (*trace_func)(N, p, k);
  }
  return res;
}

GEN trace_primes(long N, long k, long prec, int newspace, long start)
{
  GEN p_list = primes0(mkvec2(mkintn(1,start), nextprime(mkintn(1,prec))));
  long num_primes = lg(p_list);
  // adding trace of the identity (dimension of the space)
  // GEN res = cgetg(num_primes-1, t_VEC);
  GEN res = cgetg(num_primes, t_VEC);
  long p;
  GEN NK = mkvec2(mkintn(1,N),mkintn(1,k));
  GEN f = mftraceform(NK,newspace ? 0 : 1);

  gel(res, 1) = mfcoef(f, 1);
  
  for (long idx = 1; idx < num_primes - 1; idx++)
  {
    p = gtos(gel(p_list, idx));
    // gel(res, idx) = mfcoef(f, p);
    gel(res, idx+1) = mfcoef(f, p);
  }
  return res;
}

GEN traceALupto(long N, long k, long prec, long start)
{
  GEN res = cgetg(prec+1-start, t_VEC);
  long i = 1;
  
  if (start == 0) {
    // We add a 0 in the beginning to align with mfcoefs
    gel(res, i) = gen_0;
    i++;
  }
  for (; i <= prec; i++)
    gel(res, i) = traceAL(N, i+start-1, k);

  return res;
}

time_t timeTraceAL(long upTo, long from, long k, long num_traces, int only_primes, int newspace, long start)
{
  time_t start_time = time(NULL);
  long prec;
  GEN NK;
  GEN res, f, coefs, all_coefs;

  char filename[80];
  FILE* outfile;

#ifdef DEBUG
  printf("In timeTraceAL, with upTo = %ld\n", upTo);
#endif // DEBUG
  
  for (long N = from; N < upTo; N++) 
  {
    if (N == 0) continue;
    // p = gtos(gel(p_list, idx));
    // sprintf(p_str, "%d", p);
    snprintf(filename, 30, "data/traces_%ld_%ld.m", k, N);
#ifdef DEBUG
    printf("output directed to file %s\n", filename);
#endif // DEBUG
    if (num_traces == -1) {
      prec = maxuu((k*N+11) / 12, 1000);
      prec = maxuu(prec, 30*sqrt(N));
    }
    else {
      prec = num_traces;
    }
    if (only_primes)
      res = traceALprimes(N, k, prec+1, newspace, start);
    else
      res = traceALupto(N, k, prec+1, start);

    if (only_primes)
      coefs = trace_primes(N, k, prec+1, newspace, start);
    else {
      NK = mkvec2(mkintn(1,N),mkintn(1,k));
      f = mftraceform(NK,1);
      all_coefs = mfcoefs(f, prec, 1);
      coefs = cgetg(prec+1-start, t_VEC);
#ifdef DEBUG
      printf("getting trace form coefficients, start = %ld...\n", start);
#endif // DEBUG
      for (long i = start+1; i <= prec; i++) {
	gel(coefs, i-start) = gel(all_coefs, i);
      }
    }
    outfile = fopen(filename, "w");
    if (outfile == NULL)
      printf("Error! Could not open file %s for writing. skipping.\n",
	     filename);
    else
      pari_fprintf(outfile, "traces_%d := %Ps;\ntracesAL_%d := %Ps;\n",
		   N, coefs, N, res);
      //pari_printf("traces := %Ps;\ntracesAL := %Ps;\n", coefs, res);
    fclose(outfile);
  }
#ifdef DEBUG
  printf("Finished.\n");
#endif // DEBUG
  return time(NULL) - start_time;
}
