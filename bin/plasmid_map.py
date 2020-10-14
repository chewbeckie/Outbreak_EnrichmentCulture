#!/usr/bin/env python3
import pysam
import sys
import gfapy
import numpy
from scipy import stats
from scipy import special

# don't trust possible plasmids bigger than this
MAX_PLASMID_LENGTH = 300000
# don't trust taxonomic classification unless the contig is at least this big
MIN_CLASSIFICATION_LENGTH = 300000
# reads must be mapped at least this confidently
MIN_MAPQ = 50
PLSDB_PVALUE = 0.1 # rough upper bound at 1 in 10 as false positives
MIN_HIC_DISTANCE = 800

# only process this many Hi-C reads
MAX_HIC_READS = 5000000

# some thresholds for result reporting
MIN_TAXON_LINKS = 300
MIN_PLASMID_TAX_LINKS = 3
MIN_PLASMID_PLASMID_LINKS = 2

# fraction of reads that are Hi-C. TOOD: should be estimated from the dataset
hic_fraction = 0.25


gfa = gfapy.Gfa.from_file(sys.argv[1])
bam = pysam.AlignmentFile(sys.argv[2], "rb" )
kraken_class = open(sys.argv[3], 'r')
plasmid_table = open(sys.argv[4], 'r')
tax_links_output = open(sys.argv[5], 'w')
plasmid_links_output = open(sys.argv[6], 'w')

# get the reference sequences and find the small circular elements
ref_seqs = {}
is_plasmid = {}
contig_read_count = {}
for seg in gfa.segments:
    is_plasmid[seg.name] = False
    ref_seqs[seg.name] = seg.sequence
    contig_read_count[seg.name] = 0

# parse the kraken classifications
classifications = {}
for line in kraken_class:
    l = line.split('\t')
    if l[0] == 'U': continue # unclassified
    if len(ref_seqs[l[1]]) < MIN_CLASSIFICATION_LENGTH: continue # too small to classify
    name = l[2].split(' ')[0]
    classifications[l[1]] = name

line1 = plasmid_table.readline()
plasmid_names = {}
PLS_NAME_COL = 6
pname_pvalues = {}
for line in plasmid_table:
    l = line.rstrip().split('\t')
    if l[2] == 'NA': continue
    cur_pvalue = float(l[2])
    if cur_pvalue < PLSDB_PVALUE: is_plasmid[l[1]] = True
    if l[1] in pname_pvalues and pname_pvalues[l[1]] < cur_pvalue: continue
    pname_pvalues[l[1]] = cur_pvalue
    # get a new name
    ptype = l[-2]
    if ptype.find(','):
        ptype = ptype.split(',')[0]
    pname = l[PLS_NAME_COL]
    if len(pname) > 0:
        name_pos = pname.find("plasmid p")
        if name_pos > 0:
            pname = pname[name_pos+8:pname.find(",", name_pos)]
            pname += ", "+l[0]
        else:
            pname = "Similar to " + l[0]
    else:
        pname = "Similar to " + l[0]
    if ptype != 'NA':
        pname = ptype + ', ' + pname
    plasmid_names[l[1]] = pname

# initialize the contact count matrix
contacts = {}
edge_count = {}
for name1 in ref_seqs:
    edge_count[name1] = 0
    for name2 in ref_seqs:
        n1 = min(name1,name2)
        n2 = max(name1,name2)
        contacts[(n1,n2)] = 0

# Find candidate plasmids / phage
# look for edges that link a contig end to itself and that are small enough
for edge in gfa.edges:
    edge_count[edge.from_segment.name] += 1
    edge_count[edge.to_segment.name] += 1

circular = set()
for edge in gfa.edges:
    if edge.from_segment.name == edge.to_segment.name and \
        edge_count[edge.from_segment.name] == 2 and \
        is_plasmid[edge.from_segment.name] and \
        len(ref_seqs[edge.from_segment.name]) < MAX_PLASMID_LENGTH:
        circular.add(edge.from_segment.name)

print("Found "+str(len(circular))+" putative plasmids")

# parse through the bam counting links among contigs
read_iter = bam.fetch(until_eof=True)
next_read = next(read_iter)
read_count = 0
hic_count = 0
try:
    while True:
        # get a pair of aligned reads
        cur_read_1 = next_read
        next_read = next(read_iter)
        if cur_read_1.reference_name is None: continue
        if cur_read_1.mapping_quality < MIN_MAPQ: continue
        if next_read.query_name == cur_read_1.query_name:
            cur_read_2 = next_read
            next_read = next(read_iter)
        else:
            cur_read_2 = None

        if cur_read_2 is not None and cur_read_2.reference_name is not None:
            if cur_read_2.mapping_quality < MIN_MAPQ: continue
            tig1 = min(cur_read_1.reference_name,cur_read_2.reference_name)
            tig2 = max(cur_read_1.reference_name,cur_read_2.reference_name)
            if tig1 == tig2 and abs(cur_read_1.reference_start - cur_read_2.reference_start) > MIN_HIC_DISTANCE:
                hic_count += 1
            contig_read_count[cur_read_1.reference_name] += 1
            contig_read_count[cur_read_2.reference_name] += 1
            contacts[(tig1,tig2)] += 1
            read_count += 1
        if read_count > MAX_HIC_READS: break
except:
    pass
print("Processed "+str(read_count)+" reads, of which "+str(hic_count)+" are putative within-contig Hi-C reads")
hic_fraction = hic_count / read_count

# transform contig contacts to phage : taxon counts
tax_links = {}
plasmid_links = {}
tax_counts = {} # number of contacts with a taxon
plasmid_tax_counts = {} # number of contacts with a taxon
plasmid_counts = {}
for pair in contacts:
    if contacts[pair] == 0: continue
    circular_tig = None
    chromosome_tig = None
    if pair[0] in circular:
        circular_tig = pair[0]
    elif pair[0] in classifications:
        chromosome_tig = pair[0]
    if pair[1] in circular:
        if circular_tig is not None:
            # two circular contigs with contact -- build a plasmid network
            if pair not in plasmid_links: plasmid_links[pair]=0
            plasmid_links[pair] += contacts[pair]
            if not pair[0] in plasmid_counts: plasmid_counts[pair[0]] = 0
            plasmid_counts[pair[0]] += contacts[pair]
            if not pair[1] in plasmid_counts: plasmid_counts[pair[1]] = 0
            plasmid_counts[pair[1]] += contacts[pair]
            continue
        else:
            circular_tig = pair[1]
    elif pair[1] in classifications:
        if chromosome_tig is not None:
            if not classifications[pair[0]] in tax_counts: tax_counts[classifications[pair[0]]] = 0
            tax_counts[classifications[pair[0]]] += contacts[pair]
            if not classifications[pair[1]] in tax_counts: tax_counts[classifications[pair[1]]] = 0
            tax_counts[classifications[pair[1]]] += contacts[pair]
            continue
        chromosome_tig = pair[1]

    if circular_tig is None or chromosome_tig is None:
        continue # not a valid link

    taxpair = (circular_tig,classifications[chromosome_tig])
    if not taxpair in tax_links: tax_links[taxpair] = 0
    tax_links[taxpair] += contacts[pair]
    if not taxpair[1] in tax_counts: tax_counts[taxpair[1]] = 0
    tax_counts[taxpair[1]] += contacts[pair]
    if not taxpair[0] in plasmid_tax_counts: plasmid_tax_counts[taxpair[0]] = 0
    plasmid_tax_counts[taxpair[0]] += contacts[pair]
    if not taxpair[0] in plasmid_counts: plasmid_counts[taxpair[0]] = 0
    plasmid_counts[taxpair[0]] += contacts[pair]

print(tax_counts)
print(read_count)

used_plasmid = {}
used_tax = {}
used_taxplas = {}
for link in tax_links:
    if tax_counts[link[1]] < MIN_TAXON_LINKS: continue # ignore this one
    if plasmid_tax_counts[link[0]] < MIN_PLASMID_TAX_LINKS: continue # ignore
    if link[1] == 'Enterobacterales': continue
    if link[1] == 'Enterobacteriaceae': continue
    if link[1] == 'Proteobacteria': continue
    if tax_links[link] < 2: continue # don't believe anything seen only once
    # compute significance
    expected = (hic_fraction * plasmid_counts[link[0]] * tax_counts[link[1]]) / (2*read_count)
    log_cdf = stats.poisson.logcdf(expected,tax_links[link])
    tax_links_output.write(plasmid_names[link[0]]+'\t'+link[1])
    tax_links_output.write('\t'+str(-log_cdf)+'\t'+str(tax_links[link])+'\t'+str(expected)+'\n')
    used_plasmid[link[0]]=1
    used_tax[link[1]]=1
    used_taxplas[link]=1

for plasmid in used_plasmid:
    for taxon in used_tax:
        if not (plasmid,taxon) in used_taxplas:
            tax_links_output.write(plasmid_names[plasmid]+'\t'+taxon+'\t0\t0\t0\n')

for link in plasmid_links:
    if plasmid_links[link] < MIN_PLASMID_PLASMID_LINKS: continue
    # TODO: compute link significance under a Poisson model
    plasmid_links_output.write('\t'.join(link)+'\t'+str(plasmid_links[link])+'\n')
