# ercc1-mouse

- First QC: ```sh batch_qc.sh -i /data/public/jkoubele/ercc1/FASTQ -o /data/public/jkoubele/ercc1/QC_before_trimming```
- Detecting adapters: ```sh batch_detect_adapters.sh -i /data/public/jkoubele/ercc1/FASTQ -o /data/public/jkoubele/ercc1/detected_adapters```
- Alignment: ```sh batch_align.sh -i /data/public/jkoubele/ercc1/FASTQ -o /data/public/jkoubele/ercc1/STAR_output -g /data/public/jkoubele/reference_genomes/GRCm39 -b /data/public/jkoubele/ercc1/barcodes```  