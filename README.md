# BOVERI-561
Extract indel calls from commercial samples.

The list of commercial samples run IDs is in BOVERI-532. The corresponding input
file is in `data/BOVERI-532.csv`.

## Run-specific results

The code in this repo generates one TSV file per run, containing the list of
all (almost) unfiltered indels. The only filters that have been applied are
- discarding alignments with at least 5 gaps from variants support,
- filtering out variants supported by an alignment in a low-quality region
  defined as a region of 3 bases pre and post breakpoint with an average quality
  within a read cluster that is at most 20 for all bases of the considered
  region.

Each entry in the TSV file for a run contains the usual fields of a variant
(chromosome, position, reference sequence, alternate sequence, VAF) plus the
following fields
- sample: sample ID where the variant was observed
- run_id: run ID containing this sample
- run_name: run name
- v_type: DEL, INS, DELINS, MNV
- score: penalty associated to the variant (high penalty implies low confidence)
- complexity: penalty associated to the complexity of the sequence around the
              breakpoint (5 bases pre and post)
- support: support penalty associated to the size of the largest cluster
           supporting the variant
- overlap: overlap penalty associated to the VAF of overlapping variants
- control: closest VAF of the same variant in a control sample (0.0 if not
           present in a control sample)
- cov: coverage of the amplicon(s) by merged reads (a merged read corresponds
       to two reads)
- alt_cov: coverage of the variant by merged reads
- max_cov: size (in number of merged reads) of the largest supporting cluster
- repeats: WT1:WT2:WT3:WT4,V1:V2:V3:V4
    WT1 = length of the repeat unit of the reference sequence
    WT2 = number of copies of the repeat unit in the reference sequence
    WT3 = number of copies of the repeat unit left of variant breakpoint
    WT4 = number of copies of repeat unit right of variant breakpoint
    V1, V2, V3,  V4: same for alternate sequence
- annotation: snpEff annotation

## Aggregated results

The directory `results` contains a list of files prefixed by
`v51NextSeq_commercial` and `v51NextSeq_clinical` obtained by the script
`BOVERI-568.sh` that shows the result on series of clinical and
commercial samples dilutions.
