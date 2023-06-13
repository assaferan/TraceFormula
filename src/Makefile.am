bin_PROGRAMS = traceAL_sta traceALbatch_sta traceAL_dyn traceALbatch_dyn

ALL_CFLAGS = -I. -I${PARI_PATH}/include -Wall -Werror -ffp-contract=off -fno-strict-aliasing

OPT_CFLAGS = -O3 -DNDEBUG -msse3 -mavx
# DBG_CFLAGS = -g -DDEBUG -fsanitize=address
OPT_CFLAGS += -fomit-frame-pointer

if LINUX
	OPT_CFLAGS += -fvar-tracking-assignments-toggle
	OPT_CFLAGS += -funroll-all-loops
        # DBG_CFLAGS += -pg
endif

if OSX
	ALL_CFLAGS += -mmacosx-version-min=10.10
endif

OPT_CFLAGS += $(ALL_CFLAGS)
# DBG_CFLAGS += ${ALL_CFLAGS}

# DBG_FULL_CFLAGS = -DDEBUG_LEVEL_FULL
# DBG_FULL_CFLAGS += ${DBG_CFLAGS}

ALLLDFLAGS = -L${PARI_PATH}/lib -lpari -ffp-contract=off -fno-strict-aliasing    -Wl,--search_paths_first 

OPTLDFLAGS = ${ALLLDFLAGS}
OPTLDFLAGS += -O3

# DBG_LD_FLAGS = ${ALL_LDFLAGS}
# DBG_LD_FLAGS += -fsanitize=leak

DLCFLAGS   = -fPIC
# MODLDFLAGS = -bundle -undefined dynamic_lookup ${DLCFLAGS}
MODLDFLAGS = ${DLCFLAGS}

traceAL_sta_SOURCES = traceAL.c
traceALbatch_sta_SOURCES = traceALbatch.c

traceAL_dyn_SOURCES = traceAL.c
traceALbatch_dyn_SOURCES = traceALbatch.c

traceAL_sta_CFLAGS = ${OPT_CFLAGS}
traceAL_sta_LDFLAGS = ${OPTLDFLAGS}
traceALbatch_sta_CFLAGS = ${OPT_CFLAGS}
traceALbatch_sta_LDFLAGS = ${OPTLDFLAGS}

traceAL_dyn_CFLAGS = ${OPT_CFLAGS} ${DLCFLAGS}
traceAL_dyn_LDFLAGS = ${OPTLDFLAGS} ${MODLDFLAGS}
traceALbatch_dyn_CFLAGS = ${OPT_CFLAGS} ${DLCFLAGS}
traceALbatch_dyn_LDFLAGS = ${OPTLDFLAGS} ${MODLDFLAGS}

