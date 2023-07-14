bin_PROGRAMS = traceAL_sta traceALbatch_sta traceAL_dyn traceALbatch_dyn traceAL_sta_dbg traceALbatch_sta_dbg

ALL_CFLAGS = -I. -I${PARI_PATH}/include -Wall -Werror -ffp-contract=off -fno-strict-aliasing

OPT_CFLAGS = -O3 -DNDEBUG -msse3 -mavx
DBG_CFLAGS = -g -DDEBUG -fsanitize=address
OPT_CFLAGS += -fomit-frame-pointer

if LINUX
	OPT_CFLAGS += -fvar-tracking-assignments-toggle
	OPT_CFLAGS += -funroll-all-loops
        DBG_CFLAGS += -pg
endif

if OSX
	ALL_CFLAGS += -mmacosx-version-min=10.10
endif

OPT_CFLAGS += $(ALL_CFLAGS)
DBG_CFLAGS += ${ALL_CFLAGS}

DBG_FULL_CFLAGS = -DDEBUG_LEVEL_FULL
DBG_FULL_CFLAGS += ${DBG_CFLAGS}

ALL_LD_FLAGS = -L${PARI_PATH}/lib -lpari -ffp-contract=off -fno-strict-aliasing -lm -DSO

OPT_LD_FLAGS = ${ALL_LD_FLAGS}
OPT_LD_FLAGS += -O3

DBG_LD_FLAGS = ${ALL_LD_FLAGS}
# DBG_LD_FLAGS += -fsanitize=leak

DLCFLAGS   = -fPIC
# MODLDFLAGS = -bundle -undefined dynamic_lookup ${DLCFLAGS}
MODLDFLAGS = ${DLCFLAGS}

traceAL_sta_SOURCES = traceAL.c traceAL_tools.c
traceALbatch_sta_SOURCES = traceALbatch.c traceAL_tools.c

traceAL_sta_dbg_SOURCES = traceAL.c traceAL_tools.c
traceALbatch_sta_dbg_SOURCES = traceALbatch.c traceAL_tools.c

traceAL_dyn_SOURCES = traceAL.c traceAL_tools.c
traceALbatch_dyn_SOURCES = traceALbatch.c traceAL_tools.c

traceAL_sta_CFLAGS = ${OPT_CFLAGS}
traceAL_sta_LDFLAGS = ${OPT_LD_FLAGS}
traceALbatch_sta_CFLAGS = ${OPT_CFLAGS}
traceALbatch_sta_LDFLAGS = ${OPT_LD_FLAGS}

traceAL_sta_dbg_CFLAGS = ${DBG_CFLAGS}
traceAL_sta_dbg_LDFLAGS = ${DBG_LD_FLAGS}
traceALbatch_sta_dbg_CFLAGS = ${DBG_CFLAGS}
traceALbatch_sta_dbg_LDFLAGS = ${DBG_LD_FLAGS}

traceAL_dyn_CFLAGS = ${OPT_CFLAGS} ${DLCFLAGS}
traceAL_dyn_LDFLAGS = ${OPT_LD_FLAGS} ${MODLDFLAGS}
traceALbatch_dyn_CFLAGS =  ${DLCFLAGS} ${OPT_CFLAGS} # ${DBG_FULL_CFLAGS}
traceALbatch_dyn_LDFLAGS =  ${MODLDFLAGS} ${OPT_LD_FLAGS} # ${DBG_LD_FLAGS}

