import sys, os
from math import log
from collections import defaultdict

import leafify, asm_node

# TODO: attempt to find these or error out
SCRIPT_LOCATION = os.path.abspath(os.path.dirname(__file__))
MINIMAP = SCRIPT_LOCATION + '/minimap2/minimap2'
MINIASM = SCRIPT_LOCATION + '/miniasm/miniasm'
RACON = SCRIPT_LOCATION  + '/racon/build/bin/racon'

if len(sys.argv) != 4 and len(sys.argv) != 5:
    print("Usage: %s <reads> <alignments> <min leaf size> [alignment bin size]" % sys.argv[0], file=sys.stderr)
    sys.exit(1)

READS = sys.argv[1]
ALIGNMENTS = sys.argv[2]
MIN_LEAF_SIZE = int(sys.argv[3])
BIN_SIZE = 10000
THREADS = 4

if len(sys.argv) == 5:
    BIN_SIZE = int(sys.argv[4])

topology = leafify.leafify(READS, ALIGNMENTS, MIN_LEAF_SIZE, BIN_SIZE)

for r, nleaves in enumerate(topology):
    # Assemble first level (minimap-miniasm-racon)
    for i in range(nleaves):
        prefix = 'tmp.%d.%d.%d' % (r, i, 0)

        # TODO: check for errors
        os.system(MINIMAP + ' -t %d -x ava-pb %s %s > %s' % (THREADS, prefix + '.fa', prefix + '.fa', prefix + '.paf'))

        os.system(MINIASM + ' -f %s %s > %s' % (prefix + '.fa', prefix + '.paf', prefix + '.gfa'))
        os.system('awk \'/^S/{print ">%d.%d.0."$2"\\n"$3}\' %s > %s' % (r, i, prefix + '.gfa', prefix + '.asm.unpolished.fa'))

        os.system(MINIMAP + ' -t %d -x map-pb %s %s > %s' % (THREADS, prefix + '.asm.unpolished.fa', prefix + '.fa', prefix + '.map.paf'))
        os.system(RACON + ' -t %d %s %s %s > %s' % (THREADS, prefix + '.fa', prefix + '.map.paf', prefix + '.asm.unpolished.fa', prefix + '.asm.fa'))

    # Assemble the rest of the hierarchy (greedy)
    j = nleaves
    for l in range(1, int(log(nleaves, 2))+2):
        k = 0
        for i in range(0, j, 2):
            prefix = 'tmp.%d.%d.%d' % (r, k, l)

            os.system('cat %s %s > %s' % ('tmp.%d.%d.%d.asm.fa' % (r, i, l-1), 'tmp.%d.%d.%d.asm.fa' % (r, i+1, l-1), prefix + '.fa'))
            os.system(MINIMAP + ' -t %d -x ava-pb %s %s > %s' % (THREADS, prefix + '.fa', prefix + '.fa', prefix + '.paf'))
            asm_node.assemble(prefix + '.paf', prefix + '.fa', prefix + '.asm.fa')
            
            k += 1
        j = k

assemblies = ['tmp.%d.0.%d.asm.fa' % (r, int(log(nleaves, 2))+1) for r, nleaves in enumerate(topology)]
os.system('cat' + ' '.join(assemblies) + ' > final.asm.fa')
