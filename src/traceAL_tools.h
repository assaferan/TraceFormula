#include <pari/pari.h>
#include <stdint.h>
#include <time.h>

typedef uint64_t W64;
typedef int64_t Z64;

#ifdef __LINUX__
#else // MACOS
#endif // __LINUX__

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

GEN getcache(void);

void pari_close_mf(void);

/* D > 0; 6 * hclassno(D) (6*Hurwitz). Beware, cached value for D (=0,3 mod 4)
 * is stored at D>>1 */
ulong hclassno6w64(W64 D);

// returns 12*H
long H12(GEN D);

GEN mksintn(long l, long x);

GEN polyGegenbauer(long k, GEN t, GEN m);

long alpha(ulong n);

// In what follows we assume k >=2 is even

GEN traceAL(long N, long n, long k);

GEN traceALprimes(long N, long k, long prec, int newspace, long start);

GEN trace_primes(long N, long k, long prec, int newspace, long start);

GEN traceALupto(long N, long k, long prec, long start);

time_t timeTraceAL(long upTo, long from, long k, long num_traces, int only_primes, int newspace, long start);
