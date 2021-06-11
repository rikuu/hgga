from math import inf
import sys

# Start/end of a read
node_in = lambda n: n+'-'
node_out = lambda n: n+'+'

class OverlapGraph:
    def __init__(self):
        self.nodes = []
        self.edges = {}

    # Returns the graph and a list of contained reads
    # Parameters:
    # - PAF file
    # - minimum overlap length to consider
    # - maximum overhang to ignore in an end which overlaps with the other read
    def parse_paf(filename, min_overlap, max_overhang):
        graph = OverlapGraph()
        contained = []

        # Step 1: Find contained reads
        with open(filename, 'r') as f:
            for line in f:
                data = line.strip().split('\t')
                u, ul, us, ue = data[0], int(data[1]), int(data[2]), int(data[3])
                v, vl, vs, ve = data[5], int(data[6]), int(data[7]), int(data[8])
                strand = data[4]

                # Skip self-loops and low overlaps
                if u == v or max(ue-us, ve-vs) < min_overlap:
                    continue

                # If u is shorter than v and both overhangs in u are short,
                # then u is a contained in v and vice versa
                if ul < vl and us < max_overhang and ul-ue < max_overhang:
                    if not u in contained:
                        contained.append(u)
                elif vl <= ul and vs < max_overhang and vl-ve < max_overhang:
                    if not v in contained:
                        contained.append(v)

        # Step 2: Construct graph leaving out contained reads
        with open(filename, 'r') as f:
            for line in f:
                data = line.strip().split('\t')
                u, ul, us, ue = data[0], int(data[1]), int(data[2]), int(data[3])
                v, vl, vs, ve = data[5], int(data[6]), int(data[7]), int(data[8])
                strand = data[4]

                # Skip self-loops and low overlaps
                if u == v or max(ue-us, ve-vs) < min_overlap:
                    continue

                # Skip contained reads
                if u in contained or v in contained:
                    continue

                # Condition for the ovelap to extend (i.e. be useful)
                # That is: the reads don't have a long overhang in the same end
                if strand == '+' and \
                        ((us > max_overhang and vs > max_overhang) or
                        (ul-ue > max_overhang and vl-ve > max_overhang)):
                    continue

                if strand == '-' and \
                        ((us > max_overhang and vl-ve > max_overhang) or
                        (ul-ue > max_overhang and vs > max_overhang)):
                    continue

                # Add u to graph if it is not there already
                if not u in graph.nodes:
                    graph.nodes.append(u)
                
                    graph.edges[node_in(u)] = []
                    graph.edges[node_out(u)] = []

                # Add u to graph if it is not there already
                if not v in graph.nodes:
                    graph.nodes.append(v)
                
                    graph.edges[node_in(v)] = []
                    graph.edges[node_out(v)] = []

                # If overhang in the end is longer than in the start,
                # then the start of the read is involved in the overlap (and vice versa)
                up = node_in(u) if ul - ue > us else node_out(u)
                vp = node_in(v) if vl - ve > vs else node_out(v)

                # Add the edges
                graph.edges[up].append((vp, us, ue, ul, vs, ve, vl, strand))
                graph.edges[vp].append((up, vs, ve, vl, us, ue, ul, strand))

            print('graph: %i nodes, %i edges' % (len(graph.nodes), len(graph.edges)), file=sys.stderr)
            # print('contained: [' + ''.join([str(elem) for elem in contained]) + ']', file=sys.stderr)

        return graph, contained

    # Find a maximum path starting from node start.
    # The path will always start with an edge between in and out nodes of the same read.
    # After that it will alternate between overlap edges and edges between in and out nodes of the same read.
    # Parameters:
    # - the node to start from
    # - the set of visited nodes which must not be included in the path
    # Returns: the maximum path and the updated list of visited nodes
    def max_path(self, start, visited):
        path = []
        n = start
        
        # move from n to the other node representing the same read
        nout = n[:-1] + '+' if n[-1] == '-' else n[:-1] + '-'
        
        # Mark nodes as visited
        visited.append(n)
        visited.append(nout)

        while len(self.edges[nout]) != 0:
            # Find the possible extensions
            out = []
            for v, us, ue, ul, vs, ve, vl, strand in self.edges[nout]:
                if v in visited:
                    continue

                # (node, overlap, (edge))
                out.append((v, max(ue-us, ve-vs), (nout, us, ue, ul, v, vs, ve, vl, strand)))

            # No extensions found
            if len(out) == 0:
                break

            # Choose the extension with maximum overlap
            m = max(out, key=lambda x: x[1])
            n = m[0]
            
            # move from n to the other node representing the same read
            nout = n[:-1] + '+' if n[-1] == '-' else n[:-1] + '-'
            
            # Mark nodes as visited
            visited.append(n)
            visited.append(nout)

            # Add the chosen extension to the path
            path.append(m[2])

        return path, visited

    # Returns a list of list of edges (including "self-edges")
    def max_paths(self):
        if len(self.edges) == 0:
            return []
        
        paths, visited = [], []
        while len(visited) < len(self.edges.keys()):
            # Find node with least incoming edges
            start, min_edges = None, inf
            for node, edges in self.edges.items():
                if node in visited:
                    continue

                if len(edges) < min_edges:
                    start = node
                    min_edges = len(edges)

            path, visited = self.max_path(start, visited)
            paths.append(path)

        return paths

    # Transfroms list of edges to string operations
    @staticmethod
    def simplify_path(path):
        simplified = []
        for i in range(len(path)):
            u, us, ue, ul, v, vs, ve, vl, strand = path[i]
            if u[:-1] != v[:-1]:
                # Add starting substring
                if len(simplified) == 0:
                    simplified.append((u[:-1], 0, ul, u[-1]))

                # Check if u has an overhang in this end as well. If so, shorten u accordingly
                if u[-1] == '-':
                    if us > 0:
                        uentry = simplified[-1]
                        simplified[-1] = (uentry[0], us, uentry[2], uentry[3])
                        # print('Clipping %s bases' % us, file=sys.stderr)
                else:
                    if ul - ue > 0:
                        uentry = simplified[-1]
                        simplified[-1] = (uentry[0], uentry[1], ue, uentry[3])
                        # print('Clipping %s bases' % (ul-ue), file=sys.stderr)

                # Add the overhang of v
                if v[-1] == '-':
                    simplified.append((v[:-1], ve, vl, '+'))
                else:
                    simplified.append((v[:-1], 0, vs, '-'))

        return simplified
