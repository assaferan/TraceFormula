#! /usr/local/bin/sage -python

from glob import glob
import os
import time

k = 2
min_N = 1
max_N = 100000
batch_size = 50
num_in_command = 50
fnames = glob("data/traces_*")
covered = {int(fname.split('_')[-1].split('.')[0]) for fname in fnames}
not_covered = sorted([x for x in range(min_N, max_N) if x not in covered])
tmp = {x // num_in_command for x in not_covered}
tmp = sorted([t for t in tmp])
commands = ["./src/traceALbatch_sta " + str(min_N + num_in_command*i) + " " + str(min_N + num_in_command*(i + 1)) + " " + str(k) + " &" for i in tmp]
batches = [commands[i:i+batch_size] for i in range(0,len(commands),batch_size)]
for batch in batches:
    N_batch = int(batch[-1].split()[-3])
    levels_tiny = {N for N in range(min_N, N_batch)}
    remain_tiny = sorted([x for x in levels_tiny if x not in covered])
    # print("Running the following commands:", batch)
    tmp = [os.system(cmd) for cmd in batch]
    while len(remain_tiny) > 0:
    	  fnames = glob("data/traces_*")
    	  covered = {int(fname.split('_')[-1].split('.')[0]) for fname in fnames}
    	  remain_tiny = sorted([x for x in levels_tiny if x not in covered])
    	  print(len(levels_tiny) - len(remain_tiny), "/", len(levels_tiny))
    	  time.sleep(10)
    print("Finished up to ", N_batch)
