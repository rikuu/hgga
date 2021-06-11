import sys, pysam
from collections import defaultdict
from itertools import zip_longest
from math import ceil

def grouper(iterable, n, fillvalue=None):
    args = [iter(iterable)] * n
    return zip_longest(*args, fillvalue=fillvalue)

def parse_alignments(filename, bin_size):
    alignments = defaultdict(list)
    with pysam.AlignmentFile(filename, "r") as f:
        for read in f:
            if not read.is_unmapped:
                alignments[read.reference_name].append((read.reference_start // bin_size, read.query_name))
    
    return alignments

def parse_colors(filename):
    colors = defaultdict(list)
    with open(filename, 'r') as f:
        for line in f:
            data = line.strip().split()
            read, start, end = data[0], int(data[1]), int(data[2])

            # Color 0 denotes uncolored read
            if start == 0 or end == 0:
                continue

            # low bits is bin, high bits is contig
            sb, eb = start & 0xFFFF, end & 0xFFFF
            st, et = start >> 16, end >> 16

            # assert st == et

            # when st != et, coloring is chimeric in contig space
            # this is fine in most cases, however arithmetic on the coloring is sus
            # this could be handled better than throwing away high bits of the color

            mid = start + ((eb - sb) // 2)
            colors[st].append((mid, read))

    return colors

# Returns a list of trees and a count of leaves in each tree
def construct_hierarchies(alignments, min_leaf_size):
    # topology[i] is the number of leaves in tree i
    trees, topology = [], []
    for reads in alignments.values():
        sorting = sorted(reads, key = lambda x: x[0])

        tree = []
        leaf, leaf_size, i = 0, 0, 0
        while i < len(sorting):
            # Add next bin of reads to a leaf
            j = i
            while j < len(sorting) and sorting[j][0] == sorting[i][0]:
                tree.append((sorting[j][1], leaf))
                leaf_size += 1
                j += 1
            
            i = j
            
            # Move to next leaf
            if leaf_size >= min_leaf_size:
                leaf += 1
                leaf_size = 0

        trees.append(tree)
        topology.append(leaf + 1)

    return trees, topology

# Writes each leaf into a file
def write_leaves(trees, topology, reads_file):
    for i, tree in enumerate(trees):
        reads = defaultdict(list)
        for n, c in tree:
            reads[n].append(c)

        files = { j: open('tmp.'+str(i)+'.'+str(j)+'.0.fa', 'w') for j in range(topology[i]) }

        # NOTE: fastq -> fasta
        # TODO: support more input and output formats
        with open(reads_file) as f:
            for lines in grouper(f, 4, ''):
                name = lines[0][1:-1].split(' ')[0]
                if name in reads:
                    read = '>' + lines[0][1:] + lines[1]
                    for c in reads[name]:
                        files[c].write(read)

        for handle in files.values():
            handle.close()

def leafify(reads_file, guide_file, leaf_fraction, bin_size):
    guide = None
    if guide_file[-4:] == '.sam':
        # Using SAM alignments as guide
        print('SAM leafify, ' + guide_file, file=sys.stderr)
        guide = parse_alignments(guide_file, bin_size)
    else:
        # Using kermit colors as guide
        print('Kermit leafify, ' + guide_file, file=sys.stderr)
        guide = parse_colors(guide_file)
    
    # Compute minimum leaf size
    num_of_reads = sum([len(reads) for reads in guide.values()])
    min_leaf_size = ceil(leaf_fraction * num_of_reads)

    # NOTE: We use half size leaves to get double the leaves for overlapping
    # We should see an average leaf size of min_leaf_size after overlapping
    trees, topology = construct_hierarchies(guide, min_leaf_size // 2)

    # Force leaves to overlap by adding left and right halves of leaves
    # to left and right neighbor
    # NOTE: This assumes leaves to be linear, i.e. left and right have meaning
    new_topology = [0 for i in range(len(trees))]
    overlapping = [[] for i in range(len(trees))]
    for t, tree in enumerate(trees):
        n_blocks = topology[t]

        if n_blocks < 4:
            new_topology[t] = n_blocks
            overlapping[t] = tree
            continue

        blocks = [[] for i in range(n_blocks)]
        for n, c in tree:
            blocks[c].append(n)
        
        overlapping[t] = []

        c = 0
        for n in blocks[0] + blocks[1] + blocks[2]:
            overlapping[t].append((n, c))

        c = 1
        for i in range(1, n_blocks-3, 2):
            for n in blocks[i] + blocks[i+1] + blocks[i+2] + blocks[i+3]:
                overlapping[t].append((n, c))
            c += 1

        for n in blocks[n_blocks-3] + blocks[n_blocks-2] + blocks[n_blocks-1]:
            overlapping[t].append((n, c))
        
        new_topology[t] = c + 1

    write_leaves(overlapping, new_topology, reads_file)
    return new_topology

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: %s <reads> <alignments> <min leaf size>" % sys.argv[0], file=sys.stderr)
        sys.exit(1)

    reads_file = sys.argv[1]
    guide_file = sys.argv[2]
    min_leaf_size = float(sys.argv[3])
    
    leafify(reads_file, guide_file, min_leaf_size, 10000)
