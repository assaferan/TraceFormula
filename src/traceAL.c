#include <pari/pari.h>
#include <time.h>

#include "traceAL_tools.h"

int main()
{
  GEN /* f, NK,*/ k, N, al, upto, from, prec;
  // long prec; // , rem;
  pari_init(10000000000,2);
  printf("N = "); N = gp_read_stream(stdin);
  printf("k = "); k = gp_read_stream(stdin);
  printf("prec = "); prec = gp_read_stream(stdin);
  time_t start = time(NULL);
  al = traceALupto(gtos(N), gtos(k), gtos(prec));
  printf("single run took %ld seconds\n", time(NULL)-start);
  pari_printf("al = %Ps\n", al);
  // NK = mkvec2(N,k);
  // f = mftraceform(NK,0);
  printf("from = "); from = gp_read_stream(stdin);
  printf("upto = "); upto = gp_read_stream(stdin);
  time_t timing = timeTraceAL(gtos(upto), gtos(from), gtos(k), -1, 1);
  printf("took %ld seconds\n", timing);
  pari_close_mf();
  pari_close();
  return 0;
}
