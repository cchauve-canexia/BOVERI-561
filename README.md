# BOVERI-561
Extract indel calls from commercial samples.

The list of commercial samples run IDs is in BOVERI-532.

The code in this repo generates one TSV file per run, containing the list of all unfiltered indels.

The only filters that have been applied were
- discarding alignments with at least 5 gaps,
- filtering out variants supported by an alignment in a low-quality region defined as a region of 3 bases 
  pre and post breakpoint with a maximum quality of 20.
  
Each entry in the TSV file contains the usual fields of a variant (chromosome, position, reference sequence,
alternate sequence, VAF) plus the following fields
- sample: sample ID where the variant was observed
- run_id: run ID containing this sample
- run_name: run name
- v_type: DEL, INS, DELINS, MNV
- score: penalty associated to the variant (high penalty implies low confidence)
- complexity: penalty associated to the complexity of the sequence around the breakpoint (5 bases pre and post)
- support: support penalty associated to the size of the largest cluster supporting the variant
- overlap: overlap penalty associated to the VAF of overlapping variants
- control: closest VAF of the same variant in a control sample (0.0 if not present in a control sample)
- cov: coverage of the amplicon(s) by merged reads (a merged read corresponds to two reads)
- alt_cov: coverage of the variant by merged reads
- max_cov: size (in number of merged reads) of the largest supporting cluster
- annotation: snpEff annotation


