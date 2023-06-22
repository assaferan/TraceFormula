#include <pari/pari.h>
#include <time.h>

#include "traceAL_tools.h"

int
main(int argc, char* argv[])
{
  if ((argc < 4) || (argc > 6))
    {
      printf("Incorrect number of arguments.\n");
      printf("Usage: %s <from> <to> <weight> (<num_traces>) (<only_primes>) \n", argv[0]);
      printf("computes traces of S_k(N)^+ for level N between from and to, weight k and for Hecke operators up to the greatest between the Sturm bound, 30 sqrt(p) and 1000 if not specified.\n");
      return -1;
    }
  long from = atoi(argv[1]);
  long upto = atoi(argv[2]);
  long k = atoi(argv[3]);
  long num_traces = -1;
  int only_primes = 0;
  if (argc >= 5)
    num_traces = atoi(argv[4]);
  if (argc == 6)
    only_primes = atoi(argv[5]);

  pari_init(10000000000,2);
  timeTraceAL(upto, from, k, num_traces, only_primes);
  // time_t timing = timeTraceAL(upto, from, k);
  // printf("took %ld seconds\n", timing);
  pari_close_mf();
  pari_close();
  return 0;
}
