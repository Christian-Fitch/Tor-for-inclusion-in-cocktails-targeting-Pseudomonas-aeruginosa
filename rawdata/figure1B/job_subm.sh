#!/bin/bash
#SBATCH --job-name=vc3_cf
#SBATCH --export=ALL
#SBATCH -D .
#SBATCH --partition=pq # submit to the parallel queue
#SBATCH --time=96:00:00 # maximum walltime for the job
#SBATCH --account=Research_Project-172179 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --output=vc3_cf.log
source activate vcontact3

vcontact3 run --nucleotide cat_all_genomes.fasta --output vc3_fullprokdb_cf --db-domain "prokaryotes" --exports cytoscape profiles completeness --db-path ~/db/
