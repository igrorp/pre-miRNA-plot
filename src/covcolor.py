import HTSeq
import collections

gfffile = '/home/igror/Documents/lgpp/premirnaplot/premirnagit/src/covcolor/R265_2020-01-30.gff3'
bamfile = 'src/covcolor/WT_YNB1_namesorted.bam'

# for line in open(gfffile).readlines():

#     try:

#         chr, init, startalignfile = []
# for alig in HTSeq.BAM_Reader(bamfile):
    
#     if alig.aligned:

#         alignfile.append(alig)

# cvg = HTSeq.GenomicArray("auto", stranded=True, typecode='i')
    
#     except ValueError:
#         print('Error! Could not parse GFF file')
#         quit()

# alignfile = []
# for alig in HTSeq.BAM_Reader(bamfile):
    
#     if alig.aligned:

#         alignfile.append(alig)

# cvg = HTSeq.GenomicArray("auto", stranded=True, typecode='i')

# for alignment in aligfile:

#     #print(alignment)
    
#     if alignment.aligned:
    
#         cvg[alignment.iv] += 1

gff_file = HTSeq.GFF_Reader(gfffile)
features = HTSeq.GenomicArrayOfSets( "auto", stranded=True )

for feature in gff_file:
    features[feature.iv] += feature.attr["Name"]


almnt_file = HTSeq.BAM_Reader(bamfile)
counts = collections.Counter( )

for bundle in HTSeq.pair_SAM_alignments( almnt_file, bundle=True ):

   if len(bundle) != 1:
      continue  # Skip multiple alignments

   first_almnt, second_almnt = bundle[0]  # extract pair

   if not first_almnt.aligned and second_almnt.aligned:
      counts["_unmapped"] += 1
      continue
   gene_ids = set()
   for iv, val in features[left_almnt.iv].steps():
      gene_ids |= val
   for iv, val in features[right_almnt.iv].steps():
      gene_ids |= val
   if len(gene_ids) == 1:
      gene_id = list(gene_ids)[0]
      counts[ gene_id ] += 1
   elif len(gene_ids) == 0:
      counts[ "_no_feature" ] += 1
   else:
      counts[ "_ambiguous" ] += 1

for gene_id in counts:
   print(gene_id, counts[gene_id])
