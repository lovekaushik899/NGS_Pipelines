# NGS Analysis Pipeline for DNA-Seq Data

## Overview
This pipeline automates the processing of **Next-Generation Sequencing (NGS) DNA-Seq data**, starting from raw FASTQ files and ending with sorted and indexed BAM files. It follows a structured approach with **separate Conda environments** for each tool to avoid dependency conflicts.  

---

## Pipeline Workflow

1. **Quality Control of Raw Reads**  
   - Uses **FastQC** to analyze the quality of raw reads.  
   - MultiQC aggregates reports from FastQC.

2. **Trimming of Low-Quality Reads & Adapters**  
   - **Trim Galore** removes adapters and low-quality bases.

3. **Quality Control of Trimmed Reads**  
   - Runs another round of **FastQC & MultiQC** to ensure trimming effectiveness.

4. **Alignment of Reads to Reference Genome**  
   - Uses **BWA** and **HiSAT2** to align trimmed reads to a reference genome.

5. **Conversion of SAM to BAM Format**  
   - **Samtools** converts SAM files into BAM format.

6. **Sorting & Indexing BAM Files**  
   - **Samtools** sorts and indexes BAM files for downstream analysis.

---

## Installation & Dependencies 
 
|-----------------------------------------------------------------------------------------------------------------------|
| Tool        | Purpose                                         | Installation Command                        		      |
|-------------|-------------------------------------------------|-------------------------------------------------------|
| Conda       | Manages environments and dependencies           | Already installed in base environment      		        |
| FastQC      | Quality check of raw and trimmed reads          | `conda create -n fastqc_env -c bioconda fastqc` 	    |
| MultiQC     | Aggregation of FastQC reports                   | `conda create -n multiqc_env -c bioconda multiqc` 	  |
| Cutadapt    | Tool and a requisite for Trim Galore		| `conda create -n multiqc_env -c bioconda cutadapt`	          |
| Trim Galore | Trimming of adapters & low-quality reads       	| `conda create -n trim_env -c bioconda trim-galore` 	  |
| BWA         | Aligns reads to the reference genome (RNA+DNA)  | `conda create -n bwa_env -c bioconda bwa` 		        |
| HiSAT2      | Aligns read to reference genome (RNA-Seq)       | `conda create -n bwa_env -c hisat2`                   |
| Samtools    | Converts, sorts, and indexes BAM files          | `conda create -n samtools_env -c bioconda samtools` 	|
|-----------------------------------------------------------------------------------------------------------------------|

### Note
- Each tool is installed in **its own Conda environment** to prevent dependency conflicts.  
- Activate the required environment before running the corresponding tool:
  ```bash
  conda activate fastqc_env
  fastqc input.fastq
  conda deactivate
  ```

---

## How to Run the Pipeline
Ensure that:
1. **Raw FASTQ files** are in the working directory.  
2. The **reference genome is pre-indexed** using BWA.  
3. The sample names are listed in `reads.txt`.

### Execution Steps
```bash
bash DNA_Seq_Pipeline.sh
```

### Expected Output
- **Quality Reports** → `qc_results/`
- **Trimmed Reads** → `trimmed_reads/`
- **Aligned Reads (SAM/BAM)** → `alignment_results/`
- **Indexed BAM files** → `alignment_results/*.bai`
- **Logs** → `pipeline.log`

---

## Additional Notes
- The pipeline assumes the **genome has been indexed** beforehand.
- For genome indexing with BWA:
  ```bash
  bwa index reference_genome.fa
  ```
- To check BAM file statistics:
  ```bash
  samtools flagstat aligned_results/sample.bam
  ```

---

## Contributors
- **[Love Kaushik]**
- **Open for contributions**—feel free to fork, modify, and submit pull requests!

---

## License
This project is licensed under the **MIT License**. See [LICENSE](LICENSE) for details.
