#! /usr/bin/env python

###############################################################################
# Plot the evolution of the timings
###############################################################################

import csv
import getopt
import numpy as np
import pylab
import sys

# Read the command line arguments
argv = sys.argv[1:]
opts, args = getopt.getopt(argv, 'hi:', ['input_file='])

# Parse the arguments
for opt, arg in opts:
    if opt == '-h':
        print('python performance_plot.py -i input_file')
        sys.exit()
    elif opt in ('-i', '--input_file'):
        input_file = arg

# Extract the name of the benchmarks. We will use them as the title of our
# plots
titles = []
with open(input_file, 'r') as f:
    n_benchmarks = int(f.readline())
    for i in range(n_benchmarks):
        titles.append(f.readline())
    timing_types = f.readline().split()

# Extract the commit hashes and the timings
data = [[] for i in range(n_benchmarks)]
commit_hash = []
build_number = []
with open(input_file, 'r') as f:
    reader = csv.reader(f)
    for i in range(n_benchmarks + 2):
        reader.next()
    counter = 0
    for row in reader:
        if counter == 0:
            commit_hash.append(row)
            row = reader.next()
            build_number.append(row)
            row = reader.next()
        data[counter].append(row)
        counter = (counter + 1) % n_benchmarks

# Plot the timings
for i in range(n_benchmarks):
    for j in range(len(timing_types)):
        timing = []
        for d in data[i]:
            timing.append(d[j])
        pylab.plot(timing)
        pylab.xticks(np.arange(len(build_number)), build_number, rotation=90)
        pylab.ylabel('Time ($\mu$s)')
        pylab.title(titles[i])
        pylab.grid(True, which="both")
        pylab.tight_layout()
        pylab.savefig(timing_types[j] + '_' + str(i) + '.png')
        pylab.clf()
        pylab.cla()
