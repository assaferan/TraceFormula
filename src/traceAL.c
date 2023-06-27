#include <pari/pari.h>
#include <time.h>

#include "traceAL_tools.h"

int main(int argc, char* argv[])
{
  GEN al;
  long k, N, upto, from, prec, start;
  int only_primes, newspace;

  if ((argc < 5) || (argc > 9))
    {
      printf("Incorrect number of arguments.\n");
      printf("Usage: %s <level> <weight> <num_traces> <new> (<only_primes>) (<start>) (<from> <upto>) \n", argv[0]);
      printf("computes traces of S_k(N)^+ for level N between from and to, weight k and for Hecke operators up to the greatest between the Sturm bound, 30 sqrt(p) and 1000 if not specified.\n");
      return -1;
    }
  N = atoi(argv[1]);
  k = atoi(argv[2]);
  prec = atoi(argv[3]);
  newspace = atoi(argv[4]);
  newspace = (newspace ? 1 : 0); // making sure it is either 1 or 0
  if (newspace)
    only_primes = 1;
  else
    only_primes = atoi(argv[5]);
  start = 0;
  if (argc >= 6) {
    start = atoi(argv[5 + (1 - newspace)]);
  }
  if (argc >= 7 + (1-newspace)) {
    from = atoi(argv[6 + (1 - newspace)]);
    upto = atoi(argv[7 + (1 - newspace)]);
  }
  
  pari_init(10000000000,2);
  
  time_t start_time = time(NULL);
  if (only_primes)
    al = traceALprimes(N, k, prec, newspace, start);
  else
    al = traceALupto(N, k, prec, start);
  
  printf("single run took %ld seconds\n", time(NULL)-start_time);
  pari_printf("al = %Ps\n", al);
  
  if (argc >= 7 + (1-newspace)) {
    from = atoi(argv[6 + (1 - newspace)]);
    upto = atoi(argv[7 + (1 - newspace)]);
    time_t timing = timeTraceAL(upto, from, k, -1, only_primes, newspace, start);
    printf("took %ld seconds\n", timing);
  }
  pari_close_mf();
  pari_close();
  return 0;
}
