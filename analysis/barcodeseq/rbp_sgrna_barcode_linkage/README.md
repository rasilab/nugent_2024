# sgRNA-barcode linkage for the RBP dual-sgRNA ReLiC library

Links each sgRNA pair in the RBP dual-sgRNA library to the two random 20-mer barcodes carried by the same clone.

## Download

- [182p1.csv.gz](https://github.com/rasilab/nugent_2024/releases/download/linkage-data/182p1.csv.gz) (39 MB gzip, 125 MB uncompressed)
- The barcode counting workflow in [../rbp_barcode_screens/scripts/run_analysis.smk](../rbp_barcode_screens/scripts/run_analysis.smk) reads the uncompressed file from `data/filtered_barcodes/182p1.csv`, so decompress it into that folder:

```bash
gunzip -c 182p1.csv.gz > data/filtered_barcodes/182p1.csv
```

## Contents

- 2,193,296 linkage rows, sorted by `linkage_count` in descending order
- 8,760 distinct inserts, with a mean of 250 barcode pairs per insert

| Column | Meaning |
| --- | --- |
| `insert_num` | sgRNA pair identifier, joined to genes through [annotations/insert_annotations.csv](annotations/insert_annotations.csv) |
| `barcode_num` | Row rank, also used to split barcodes into two groups for the random partition control |
| `barcode1` | Barcode in the reporter mRNA 3'UTR, read by amplicon lPNHS23 |
| `barcode2` | Barcode in the mCherry-T2A-puro mRNA 3'UTR, read by amplicon lPNHS24 |
| `linkage_count` | Reads supporting the linkage |

`barcode1` and `barcode2` sit in two different transcripts of the same clone. The screens quantify reporter mRNA level as `barcode1` counts normalized by `barcode2` counts, as described in Figure 4a and 4b of the [manuscript](../../../manuscript/manuscript.md).

## Do not open this file in Excel

Excel silently truncates at 1,048,576 rows, which keeps only the first 1,048,575 linkage rows:

- 47.8% of barcode pairs are retained, covering 76.7% of linkage reads
- The truncation falls at `linkage_count` = 19, so every pair below that count is lost
- The retained fraction varies from 14% to 78% across inserts

The loss is therefore not uniform across sgRNAs. Any analysis that compares barcode-level features between inserts should use the full file.
