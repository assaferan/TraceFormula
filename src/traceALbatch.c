#include <pari/pari.h>
#include <time.h>

/*
GP;install("traceAL", "GG&&", "traceAL", "./libtraceal.so");
*/

/* fa = factorization of -D > 0, return -D0 > 0 (where D0 is fundamental) */
static long
corediscs_fact(GEN fa)
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

/***************************************************************/
/*                 Generic cache handling                      */
/***************************************************************/
enum { cache_FACT, cache_DIV, cache_H, cache_D};
typedef struct {
  const char *name;
  GEN cache;
  ulong minself, maxself;
  void (*init)(long);
  ulong miss, maxmiss;
  long compressed;
} cache;

static void constfact(long lim);
static void constdiv(long lim);
static void consttabh(long lim);
static void constcoredisc(long lim);
static THREAD cache caches[] = {
{ "Factors",  NULL,  50000,    50000, &constfact, 0, 0, 0 },
{ "Divisors", NULL,  50000,    50000, &constdiv, 0, 0, 0 },
{ "H",        NULL, 100000, 10000000, &consttabh, 0, 0, 1 },
{ "CorediscF",NULL, 100000, 10000000, &constcoredisc, 0, 0, 0 }
};

static void
cache_reset(long id) { caches[id].miss = caches[id].maxmiss = 0; }
static void
cache_delete(long id) { guncloneNULL(caches[id].cache); }
static void
cache_set(long id, GEN S)
{
  GEN old = caches[id].cache;
  caches[id].cache = gclone(S);
  guncloneNULL(old);
}

/* handle a cache miss: store stats, possibly reset table; return value
 * if (now) cached; return NULL on failure. HACK: some caches contain an
 * ulong where the 0 value is impossible, and return it (typecast to GEN) */
static GEN
cache_get(long id, ulong D)
{
  // printf("In cache_get, with id = %ld, D = %lu\n", id, D);
  cache *S = &caches[id];
  const ulong d = S->compressed? D>>1: D;
  ulong max, l;

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
static GEN
cache_report(long id)
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
GEN
getcache(void)
  {
  pari_sp av = avma;
  GEN M = cgetg(6, t_MAT);
  gel(M,1) = cache_report(cache_FACT);
  gel(M,2) = cache_report(cache_DIV);
  gel(M,3) = cache_report(cache_H);
  gel(M,4) = cache_report(cache_D);
  return gerepilecopy(av, shallowtrans(M));
}

void
pari_close_mf(void)
{
  cache_delete(cache_FACT);
  cache_delete(cache_DIV);
  cache_delete(cache_H);
  cache_delete(cache_D);
}

/*************************************************************************/
/* a odd, update local cache (recycle memory) */
static GEN
update_factor_cache(long a, long lim, long *pb)
{
  const long step = 16000; /* even; don't increase this: RAM cache thrashing */
  if (a + 2*step > lim)
    *pb = lim; /* fuse last 2 chunks */
  else
    *pb = a + step;
  return vecfactoroddu_i(a, *pb);
}

/* assume lim < MAX_LONG/8 */
static void
constcoredisc(long lim)
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
      set_avma(av2); cachea = N;
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
  set_avma(av);
}

static void
constfact(long lim)
{
  pari_sp av;
  GEN VFACT = caches[cache_FACT].cache;
  long LIM = VFACT? lg(VFACT)-1: 4;
  if (lim <= 0) lim = 5;
  if (lim <= LIM) return;
  cache_reset(cache_FACT); av = avma;
  cache_set(cache_FACT, vecfactoru_i(1,lim)); set_avma(av);
}

static void
constdiv(long lim)
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
  cache_set(cache_DIV, VDIV); set_avma(av);
}

/* n > 1, D = divisors(n); sets L = 2*lambda(n), S = sigma(n) */
static void
lamsig(GEN D, long *pL, long *pS)
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
  set_avma(av); *pL = L; *pS = S;
}

/* table of 6 * Hurwitz class numbers D <= lim */
static void
consttabh(long lim)
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
        set_avma(av2); cachea = N;
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
  cache_set(cache_H, VHDH0); set_avma(av);
}

/*************************************************************************/
/* Core functions using factorizations, divisors of class numbers caches */
static GEN
myfactoru(long N)
{
  GEN z = cache_get(cache_FACT, N);
  return z? gcopy(z): factoru(N);
}

/* write n = mf^2. Return m, set f. */
static ulong
mycore(ulong n, long *pf)
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
  *pf = f; return gc_long(av,m);
}

/* write -n = Df^2, D < 0 fundamental discriminant. Return D, set f. */
static long
mycoredisc2neg(ulong n, long *pf)
{
  ulong m, D = (ulong)cache_get(cache_D, n);
  if (D) { *pf = usqrt(n/D); return -(long)D; }
  m = mycore(n, pf);
  if ((m&3) != 3) { m <<= 2; *pf >>= 1; }
  return (long)-m;
}

/* 1+p+...+p^e, e >= 1 */
static ulong
usumpow(ulong p, long e)
{
  ulong q = 1+p;
  long i;
  for (i = 1; i < e; i++) q = p*q + 1;
  return q;
}

/* Hurwitz(D0 F^2)/ Hurwitz(D0)
 * = \sum_{f|F}  f \prod_{p|f} (1-kro(D0/p)/p)
 * = \prod_{p^e || F} (1 + (p^e-1) / (p-1) * (p-kro(D0/p))) */
static long
get_sh(long F, long D0)
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
static ulong
hclassno6u_count(ulong d)
{
  ulong a, b, b2, h = 0;
  int f = 0;

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
static long
hclassno6u_2(ulong D, long D0, long F)
{
  long h;
  // printf("In hclassno6u_2, with D = %lu, D0 = %lu, F = %lu\n", D, D0, F);
  if (F == 1) h = hclassno6u_count(D);
  else
  { /* second chance */
    h = (ulong)cache_get(cache_H, -D0);
    if (!h) h = hclassno6u_count(-D0);
    h *= get_sh(F,D0);
  }
  return h;
}

/* D > 0; 6 * hclassno(D) (6*Hurwitz). Beware, cached value for D (=0,3 mod 4)
 * is stored at D>>1 */
ulong
hclassno6u(ulong D)
{
  // printf("In hclassno6u, with D = %lu\n", D);
  ulong z = (ulong)cache_get(cache_H, D);
  // printf("From cache got z = H(D) = %lu\n", z); 
  long D0, F;
  if (z) return z;
  D0 = mycoredisc2neg(D, &F);
  // printf("Fundamental disc = D0 = %ld, F = %ld\n", D0, F);
  return hclassno6u_2(D,D0,F);
}

// returns 12*H
long H12(long D)
{
  GEN uu;
  long u;
  if (D == 0) return -1;
  if (D > 0)
    switch (D % 4) {
      case 0:
      case 3:
        return 2*hclassno6u(D);
      case 1:
      case 2:
        return 0;
    }
  long is_sq = issquareall(mkintn(1,-D), &uu);
  u = gtos(uu);
  return (is_sq ? -6*u : 0); 
}

// Returns 6 * Hurwitz(D)
static ulong
hclassno6u_i(ulong D, long D0, long F)
{
  ulong z = (ulong)cache_get(cache_H, D);
  if (z) return z;
  return hclassno6u_2(D,D0,F);
}

// In what follows we assume N is prime, k = 2

long
traceAL(long N, long n)
{
  const long nN = n*N;
  const long n4N = nN << 2;
  const long n4 = n << 2;
  // Cohen does something wiser - see if it works here
  long limt, tN;
  long ret;
  long S1 = 0;

  // printf("In traceAL, N = %ld, n = %ld\n", N, n);
  limt = usqrt(n4N) / N;
  // printf("limt = %ld\n", limt);
  GEN div_nN = divisors(mkintn(1,nN));
  // pari_printf("div_nN = %Ps\n", div_nN);
  GEN div_n = divisors(mkintn(1,n));
  long num_divs_nN = lg(div_nN);
  long num_divs_n = lg(div_n);
  long phi = gtos(eulerphi(mkintn(1,N)));
  for (tN = -limt ; tN <= limt; tN++) /* t^2 < 4Nn */
  {
    long t = tN*N;
    long t2 = t*t, D = n4N - t2;
    // printf("t = %ld, D = %ld\n", t, D);
    // S1 += hclassno6u(D);
    // printf("t = %ld, D = %ld, H12(D) = %ld\n", tN, D, H12(D));
    S1 += H12(D);
  }

  if (n4 % N == 0)
  {
    long quo = n4 / N;
    for (tN = -limt ; tN <= limt; tN++) /* t^2 < 4Nn */
    {
      long t2 = tN*tN, D = quo - t2;
      // printf("t = %ld, D = %ld, H12(D) = %ld\n", tN, D, H12(D));
      //      S1 -= hclassno6u(D);
      S1 -= H12(D);
     }  
  }
  // printf("Sum of class numbers is: %ld\n", S1);
  long S2 = 0;

  // printf("num_divs_nN = %ld\n", num_divs_nN);
  for (long idx = 1; idx < num_divs_nN; idx++)
  {
    // pari_printf("div_nN[%ld] = %Ps\n", idx, gel(div_nN,idx));
    ulong d = gtos(gel(div_nN, idx));
    ulong a = nN / d;
    if ((a+d) % N == 0)
    {
      S2 += minuu(a,d);
    }
  }

  //  printf("Sum of divisors is: %ld\n", S2);
  ret = -(S1 + 12*S2*phi / N) / 24;

  for (long idx = 1; idx < num_divs_n; idx++)
  {
    ulong d = gtos(gel(div_n, idx));
    if (ugcd(N,d) == 1)
      ret +=  n / d;
  }
  return ret;
}

GEN traceALupto(long N, long prec)
{
  GEN res = cgetg(prec+1, t_VEC);
  // We add a 0 in the beginning to align with mfcoefs
  gel(res, 1) = gen_0;
  for (long i = 2; i <= prec; i++)
  {
    long trace = traceAL(N, i-1);
    gel(res,i) = gen_0;
   
    if (trace > 0)
      gel(res,i) = gadd(gel(res,i), mkintn(1,trace));
    else
      gel(res,i) = gsub(gel(res,i), mkintn(1,-trace));
     
  }
  return res;
}

time_t timeTraceAL(long upTo, long from)
{
  time_t start = time(NULL);
  long p, prec;
  GEN k = mkintn(1,2);
  GEN NK;
  GEN res, f, coefs;
  GEN p_list = primes0(mkvec2(mkintn(1,from), mkintn(1,upTo)));
  long num_primes = lg(p_list);
  // char* base_filename = "data/traces_";
  // char* suffix = ".m";
  // char p_str[10]; // we're working with up to 5 digits, so 10 should be enough
  char filename[80];
  FILE* outfile;
  
  // printf("In timeTraceAL, with upTo = %ld\n", upTo);
  // printf("num_primes = %lu\n", num_primes);
  // pari_printf("last prime = %Ps\n", gel(p_list, num_primes-1));
  for (long idx = 1; idx < num_primes-1; idx++)
  {
    p = gtos(gel(p_list, idx));
    // sprintf(p_str, "%d", p);
    sprintf(filename, "data/traces_%ld.m", p);
    // printf("output directed to file %s\n", filename);
    prec = maxuu((p+11) / 12, 1000);
    // printf("p = %ld, prec = %ld\n", p, prec);
    res = traceALupto(p, prec+1);
    NK = mkvec2(gel(p_list, idx),k);
    f = mftraceform(NK,0);
    coefs = mfcoefs(f, prec+1, 1);
    outfile = fopen(filename, "w");
    if (outfile == NULL)
      printf("Error! Could not open file %s for writing. skipping.\n",
	     filename);
    else
      pari_fprintf(outfile, "traces_%d := %Ps;\ntracesAL_%d := %Ps;\n",
		   p, coefs, p, res);
      //pari_printf("traces := %Ps;\ntracesAL := %Ps;\n", coefs, res);
    fclose(outfile);
  }
  // printf("Finished.\n");
  return time(NULL) - start;
}

int
main(int argc, char* argv[])
{
  if (argc != 3)
    {
      printf("Incorrect number of arguments.\n");
      printf("Usage: %s <from> <to>\n", argv[0]);
      printf("computes traces for primes between from and to.\n");
      return -1;
    }
  long from = atoi(argv[1]);
  long upto = atoi(argv[2]);

  pari_init(10000000000,2);
  timeTraceAL(upto, from);
  // printf("took %ld seconds\n", timing);
  pari_close_mf();
  pari_close();
  return 0;
}
