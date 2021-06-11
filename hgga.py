import sys, os
from math import log, ceil

import leafify, asm_node

# TODO: attempt to find these or error out
SCRIPT_LOCATION = os.path.abspath(os.path.dirname(__file__))
MINIMAP = SCRIPT_LOCATION + '/minimap2/minimap2'
MINIASM = SCRIPT_LOCATION + '/miniasm/miniasm'
RACON = SCRIPT_LOCATION  + '/racon/build/bin/racon'

if len(sys.argv) < 3 or len(sys.argv) > 5:
    print("Usage: %s <reads> <alignments> [number of threads] [min leaf size]" % sys.argv[0], file=sys.stderr)
    sys.exit(1)

READS = sys.argv[1]
ALIGNMENTS = sys.argv[2]
THREADS = 1
MIN_LEAF_SIZE = 0.02

# TODO: expose all of these parameters
BIN_SIZE = 10000

MIN_OVERLAP = 10000
MAX_OVERHANG = 1000
MIN_LENGTH = 1000

FINAL_MIN_OVERLAP = 10000
FINAL_MAX_OVERHANG = 1000
FINAL_MIN_LENGTH = 1000

if len(sys.argv) == 4:
    THREADS = int(sys.argv[3])

if len(sys.argv) == 5:
    MIN_LEAF_SIZE = float(sys.argv[4])

topology = leafify.leafify(READS, ALIGNMENTS, MIN_LEAF_SIZE, BIN_SIZE)

# TODO: process-level parallellism
for r, nleaves in enumerate(topology):
    # Assemble first level (minimap-miniasm-racon)
    for i in range(nleaves):
        prefix = 'tmp.%d.%d.%d' % (r, i, 0)

        # TODO: check for errors
        # TODO: support predefines for ont and hifi
        os.system(MINIMAP + ' -t %d -x ava-pb %s %s > %s' % (THREADS, prefix + '.fa', prefix + '.fa', prefix + '.paf'))

        os.system(MINIASM + ' -f %s %s > %s' % (prefix + '.fa', prefix + '.paf', prefix + '.gfa'))
        os.system('awk \'/^S/{print ">%d.%d.0."$2"\\n"$3}\' %s > %s' % (r, i, prefix + '.gfa', prefix + '.asm.unpolished.fa'))

        os.system(MINIMAP + ' -t %d -x map-pb %s %s > %s' % (THREADS, prefix + '.asm.unpolished.fa', prefix + '.fa', prefix + '.map.paf'))
        os.system(RACON + ' -t %d %s %s %s > %s' % (THREADS, prefix + '.fa', prefix + '.map.paf', prefix + '.asm.unpolished.fa', prefix + '.asm.fa'))

    # Assemble the rest of the hierarchy (greedy)
    j = nleaves
    for l in range(1, ceil(log(nleaves, 2))+1):
        k = 0
        for i in range(0, j, 2):
            prefix = 'tmp.%d.%d.%d' % (r, k, l)

            # TODO: take partial all-vs-all mappings from leaves
            os.system('cat %s %s > %s' % ('tmp.%d.%d.%d.asm.fa' % (r, i, l-1), 'tmp.%d.%d.%d.asm.fa' % (r, i+1, l-1), prefix + '.fa'))
            os.system(MINIMAP + ' -t %d -x ava-pb %s %s > %s' % (THREADS, prefix + '.fa', prefix + '.fa', prefix + '.paf'))

            # TODO: support gfa input/output
            asm_node.assemble(prefix + '.paf', prefix + '.fa', prefix + '.asm.fa',
                MIN_OVERLAP, MAX_OVERHANG, MIN_LENGTH)
            
            k += 1
        j = k

assemblies = ['tmp.%d.0.%d.asm.fa' % (r, ceil(log(nleaves, 2))) for r, nleaves in enumerate(topology)]
os.system('cat ' + ' '.join(assemblies) + ' > final.fa')

# Handle mis-separated trees by one final self-assembly
os.system(MINIMAP + ' -t %d -x ava-pb %s %s > %s' % (THREADS, 'final.fa', 'final.fa', 'final.paf'))
asm_node.assemble('final.paf', 'final.fa', 'final.asm.fa', FINAL_MIN_OVERLAP, FINAL_MAX_OVERHANG, FINAL_MIN_LENGTH)