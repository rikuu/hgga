import sys

INDEX_BIN = 300000
MAX_OVERSHOOT = 10000

def parse_markers(filenames):
    markers = {}
    for chrom, mf in enumerate(filenames):
        with open(mf) as f:
            for line in f:
                d = line.split('\t')
                if not d[0] in markers:
                    markers[d[0]] = []
                markers[d[0]].append((int(d[1]), int(d[2]), chrom))

    print('[CC] markers parsed', file=sys.stderr)
    return markers

def index_markers(markers):
    index = {}
    for contig, ms in markers.items():
        s = sorted(ms)
        index[contig] = [[] for i in range((s[-1][0] // INDEX_BIN) + 1)]
        for p, c, ch in s:
            index[contig][p // INDEX_BIN].append((p, c, ch))

    print('[CC] markers indexed', file=sys.stderr)
    return index

def parse_alignments(filename, index):
    aligned_outside = 0
    reads = {}
    with open(filename, 'r') as f:
        for line in f:
            d = line.split()
            if not d[5] in index:
                aligned_outside += 1
                continue

            if not d[0] in reads:
                reads[d[0]] = []

            # q, ql, qs, qe, strand, t, tl, ts, te, ...
            tp = [int(d[7]), int(d[8])]
            ts, te = min(tp), max(tp)
            reads[d[0]].append((d[5], ts, te))

    print('[CC] alignments parsed, %d outside' % aligned_outside, file=sys.stderr)
    return reads

def distance(s, e, p):
    if p <= e and p >= s:
        return 0
    return min(abs(p - s), abs(p - e))

def color(reads, index):
    uncolored = 0
    colors = []
    for read, alns in reads.items():
        read_colors = []
        for aln in alns:
            t, ts, te = aln

            if not t in index:
                continue
            
            idl = [(ts - MAX_OVERSHOOT) // INDEX_BIN, (te + MAX_OVERSHOOT) // INDEX_BIN]
            ids, ide = min(idl), max(idl)

            if ids >= len(index[t]):
                continue

            if ide >= len(index[t]):
                ide = len(index[t]) - 1

            indexes = []
            for i in range(ids, ide+1):
                indexes += index[t][i]

            distances = sorted([(distance(ts, te, p), c, ch) for p, c, ch in indexes])
            
            i = 0
            while i < len(distances) and distances[i][0] - distances[0][0] < MAX_OVERSHOOT:
                i += 1

            if i == 0:
                continue
            
            # TODO: majority vote for chromosome
            c = [ch << 16 | (c+1) for d, c, ch in distances[0:i+1]]
            read_colors.append((distances[0][0], min(c), max(c)))
        
        if len(read_colors) == 0:
            uncolored += 1
            continue

        best_color = sorted(read_colors, key = lambda x: (x[0], x[2] - x[1]))[0]
        colors.append((read, best_color[1], best_color[2]))

    print('[CC] %d colors, %d uncolored' % (len(colors), uncolored), file=sys.stderr)
    return colors

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: %s <mappings> <markers> [markers...]" % sys.argv[0], file=sys.stderr)
        sys.exit(1)

    paf_filename = sys.argv[1]
    markers_filenames = sys.argv[2:]
    
    markers = parse_markers(markers_filenames)
    index = index_markers(markers)
    alignments = parse_alignments(paf_filename, index)

    colors = color(alignments, index)

    for read, c1, c2 in colors:
        print('%s\t%d\t%d' % (read, c1, c2))

