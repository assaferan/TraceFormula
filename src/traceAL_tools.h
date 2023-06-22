#include <pari/pari.h>
#include <time.h>

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

static THREAD cache caches[4];

GEN getcache(void);

void pari_close_mf(void);

/* D > 0; 6 * hclassno(D) (6*Hurwitz). Beware, cached value for D (=0,3 mod 4)
 * is stored at D>>1 */
ulong hclassno6u(ulong D);

// returns 12*H
long H12(long D);

GEN mksintn(long l, long x);

GEN polyGegenbauer(long k, long t, long m);

long alpha(ulong n);

// In what follows we assume k >=2 is even

GEN traceAL(long N, long n, long k);

GEN traceALprimes(long N, long k, long prec);

GEN trace_primes(long N, long k, long prec);

GEN traceALupto(long N, long k, long prec);

time_t timeTraceAL(long upTo, long from, long k, long num_traces, int only_primes);
