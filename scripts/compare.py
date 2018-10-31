#! /usr/bin/env python

###############################################################################
# This script parses the output of google benchmark and compare the new result
# with a reference. It also extracts the timings so that they can be plotted.
###############################################################################

import csv
import getopt
import json
import os.path
import sys

# Read the command line arguments
argv = sys.argv[1:]
opts, args = getopt.getopt(argv, 'hc:b:r:n:o:t:', ['help', 'commit_hash=',
                                                   'build_number=', 'input_ref=', 'input_new=',
                                                   'output_file=', 'tolerance='])
# Default relative tolerance is 10%
tol = 0.1
# Parse the arguments
for opt, arg in opts:
    if opt in ('-h', '--help'):
        print('python compare.py -c hash -b build -r file1 -n file2 -o file3 -t float')
        print('-c commit_hash')
        print('-b build_number')
        print('-r reference_benchmark.json')
        print('-n new_benchmark.json')
        print('-o output_file.txt')
        print('-t relative tolerance (default 0.1)')
        sys.exit()
    elif opt in ('-c', '--commit_hash'):
        commit_hash = arg
    elif opt in ('-b', '--build_number'):
        build_number = arg
    elif opt in ('-r', '--input_ref'):
        ref_benchmark = arg
    elif opt in ('-n', '--input_new'):
        new_benchmark = arg
    elif opt in ('-o', '--output_file'):
        output_file = arg
    elif opt in ('-t', '--tolerance'):
        tol = arg

# Load the reference input file
with open(ref_benchmark, 'r') as f:
    ref_data = json.load(f)
# Load the new output file
with open(new_benchmark, 'r') as f:
    new_data = json.load(f)

# Check that the reference and the new run have the same number of benchmarks
assert len(ref_data['benchmarks']) == len(new_data['benchmarks'])

n_benchmarks = len(ref_data['benchmarks'])/3

# Write header. This is done only if the file does not exist:
#   - Number of benchmarks
#   - Benchmark 1 name
#   - Benchmark 2 name
#   - ...
#   - Benchmark N name
#   - Timing types
#   - Name of the reference commit
#   - Timings for the different benchmarks
timing_types = ['real_time', 'cpu_time']
if os.path.isfile(output_file) == False:
    with open(output_file, 'w') as f:
        f.write(str(n_benchmarks) + '\n')
        for i in range(n_benchmarks):
            f.write(ref_data['benchmarks'][3*i+1]['name'] + '\n')
        for t in timing_types:
            f.write(t + ' ')
        f.write('\n')
        f.write('ref\n')
        writer = csv.writer(f)
        for i in range(n_benchmarks):
            row = []
            for time in timing_types:
                ref_time = ref_data['benchmarks'][3*i+1][time]
                row.append(ref_time)
            writer.writerow(row)

# Write the commit hash (only the first seven characters), the build number,
# and the timings for the different benchmarks
build_passed = True
failing_benchmarks = []
with open(output_file, 'a') as f:
    f.write(commit_hash[0:7] + '\n')
    f.write(build_number + '\n')
    writer = csv.writer(f)
    for i in range(n_benchmarks):
        row = []
        for time in timing_types:
            # We are only interested in the median time not the average or the
            # standard deviance. We use the median because it is more stable
            # than the average.
            new_time = new_data['benchmarks'][3*i+1][time]
            ref_time = ref_data['benchmarks'][3*i+1][time]
            row.append(new_time)
            if (new_time - ref_time) / ref_time > tol:
                failing_benchmarks.append([i, time, ref_time, new_time])
                build_passed = False
        writer.writerow(row)

if build_passed == True:
    sys.exit(0)
else:
    for failure in failing_benchmarks:
        name = ref_data['benchmarks'][3*failure[0]+1]['name']
        time_name = failure[1]
        ref_time = failure[2]
        new_time = failure[3]
        increase_pct = 100.*(failure[3] - failure[2])/failure[2]
        print("Failing benchmark \"{} ({})\": reference time = {}, new time = {} [{:+.1f}%]".format(
            name, time_name, ref_time, new_time, increase_pct))
    sys.exit(1)
