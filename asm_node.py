import sys

from overlap import OverlapGraph

REV = { 'A': 'T',  'T': 'A', 'G': 'C', 'C': 'G' }
reverse_complement = lambda s: ''.join(REV[c] for c in reversed(s.upper()))

# Parse fasta file to memory
def parse_reads(filename):
    reads = {}
    with open(filename, 'r') as f:
        name, seq = '', ''
        for line in f:
            if line[0] == '>':
                if name != '':
                    reads[name] = seq
                    seq = ''
                name = line[1:].rstrip().split()[0]
            else:
                seq += line.rstrip()
        
        if name != '':
            reads[name] = seq
    
    return reads

def construct_sequence(path, graph, reads):
    seq, contained, name = '', [], ''

    if len(path) == 0:
        return seq, contained, name

    simplified = graph.simplify_path(path)
    for u, us, ue, strand in simplified:
        if strand == '+':
            seq += reads[u][us:ue]
        else:
            seq += reverse_complement(reads[u][us:ue])

        name += '%s%s%d:%d/' % (u, strand, us, ue)
        contained.append(u)

    return seq, contained, name

def assemble(mappings, reads_file, out_file, min_overlap, max_overhang, min_length):
    graph, assembled = OverlapGraph.parse_paf(mappings, min_overlap, max_overhang)
    paths = graph.max_paths()

    reads = parse_reads(reads_file)

    with open(out_file, 'w') as f:
        for path in paths:
            seq, contained_reads, name = construct_sequence(path, graph, reads)
            assembled += contained_reads
            if len(seq) > min_length:
                f.write('>%s\n%s\n' % (name, seq))

        for n, s in reads.items():
            if n in assembled:
                continue

            if len(s) < min_length:
                continue

            # print('D: %s not contained' % n, file=sys.stderr)
            f.write('>%s\n%s\n' % (n, s))

if __name__ == "__main__":
    if len(sys.argv) < 3 or len(sys.argv) > 6:
        print("Usage: %s <mappings> <reads> [min overlap] [max overhang] [min length]" % sys.argv[0], file=sys.stderr)
        sys.exit(1)

    mappings = sys.argv[1]
    reads = sys.argv[2]

    min_overlap = 10000
    if len(sys.argv) >= 4:
        min_overlap = int(sys.argv[3])

    max_overhang = 1000
    if len(sys.argv) >= 5:
        max_overhang = int(sys.argv[4])

    min_length = 1000
    if len(sys.argv) == 6:
        min_length = int(sys.argv[5])
    
    assemble(mappings, reads, reads + '.asm.fa', min_overlap, max_overhang, min_length)
