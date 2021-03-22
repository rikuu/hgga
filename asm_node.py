import sys

from overlap import OverlapGraph

MIN_OVERLAP = 10000
MAX_OVERHANG = 1000

REV = { 'A': 'T',  'T': 'A', 'G': 'C', 'C': 'G' }
reverse_complement = lambda s: ''.join(REV[c] for c in reversed(s.upper()))

def parse_reads(filename):
    reads = {}
    with open(filename, 'r') as f:
        name, seq = '', ''
        for line in f:
            if line[0] == '>':
                if name != '':
                    reads[name] = seq
                    seq = ''
                name = line[1:].rstrip()
                name = name.split()[0]
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
    
    print(path)
    print(simplified)

    for u, us, ue, strand in simplified:
        if strand == '+':
            seq += reads[u][us:ue]
        else:
            seq += reverse_complement(reads[u][us:ue])

        name += '%s%s%d:%d/' % (u, strand, us, ue)
        contained.append(u)

    return seq, contained, name

def assemble(mappings, reads_file, out_file, min_overlap=MIN_OVERLAP, max_overhang=MAX_OVERHANG):
    graph, assembled = OverlapGraph.parse_paf(mappings, min_overlap=min_overlap, max_overhang=max_overhang)
    paths = graph.max_paths()

    reads = parse_reads(reads_file)

    with open(out_file, 'w') as f:
        for i, path in enumerate(paths):
            seq, contained_reads, name = construct_sequence(path, graph, reads)
            if seq == '':
                continue

            f.write('>%s\n%s\n' % (name, seq))
            assembled += contained_reads
            print(name)

        # TODO: Actually find contained reads, maybe?
        for n, s in reads.items():
            if n in assembled:
                continue

            print('D: %s not contained' % n)
            f.write('>%s\n%s\n' % (n, s))

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: %s <mappings> <reads>" % sys.argv[0], file=sys.stderr)
        sys.exit(1)

    mappings = sys.argv[1]
    reads = sys.argv[2]
    
    assemble(mappings, reads, reads + '.asm.fa')
