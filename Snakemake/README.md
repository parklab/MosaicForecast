## before running:
1. Install snakemake by using conda:
	source {your_dir}/miniconda3/etc/profile.d/conda.sh
	conda activete MF
	conda install -c bioconda -c conda-forge snakemake
2. Soft-link your case .bam and .bai files under "bam_links_case", and soft-link your panel-of-normal file under "bam_links_PON".
3. Change all directory names and fasta file names in the "Snakefile".

## generate your pipeline figure:
snakemake --dag -np |dot -Tpng > pipeline.png

## dry run test:
snakemake -np

## actual run:
snakemake

## run on cluster:

# on LSF:
source {your_dir}/miniconda3/etc/profile.d/conda.sh
conda activate MF
snakemake --unlock
snakemake --rerun-incomplete -j {job_num} --cluster-config LSF.json --cluster "bsub -n 4 --latency-wait 120 -W 60:00 -M 15G -q {queue} -o logs/%J.out -e logs/%J.err"

# on slurm:
source {your_dir}/miniconda3/etc/profile.d/conda.sh
conda activate MF
snakemake --unlock
snakemake --rerun-incomplete -j {job_num} --cluster-config cluster.json --cluster "sbatch -p {queue} --account={accountID} -c 1 -t 2-12:00 --latency-wait 120 --mem=5000 -o logs/%j.out -e logs/%j.err "
