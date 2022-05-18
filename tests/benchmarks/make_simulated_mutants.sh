#for i in {1..20}; 
#do
#	python mutate_seqs.py -i seqs/hg38.experiment.$i.fasta -r 0.01 -o seqs/hg38.experiment.$i.001.fasta
#	python mutate_seqs.py -i seqs/hg38.experiment.$i.fasta -r 0.02 -o seqs/hg38.experiment.$i.002.fasta
#	python mutate_seqs.py -i seqs/hg38.experiment.$i.fasta -r 0.03 -o seqs/hg38.experiment.$i.003.fasta
#	python mutate_seqs.py -i seqs/hg38.experiment.$i.fasta -r 0.04 -o seqs/hg38.experiment.$i.004.fasta
#	python mutate_seqs.py -i seqs/hg38.experiment.$i.fasta -r 0.05 -o seqs/hg38.experiment.$i.005.fasta
#	python mutate_seqs.py -i seqs/hg38.experiment.$i.fasta -r 0.06 -o seqs/hg38.experiment.$i.006.fasta
#	python mutate_seqs.py -i seqs/hg38.experiment.$i.fasta -r 0.07 -o seqs/hg38.experiment.$i.007.fasta
#	python mutate_seqs.py -i seqs/hg38.experiment.$i.fasta -r 0.08 -o seqs/hg38.experiment.$i.008.fasta
#	python mutate_seqs.py -i seqs/hg38.experiment.$i.fasta -r 0.09 -o seqs/hg38.experiment.$i.009.fasta
#	python mutate_seqs.py -i seqs/hg38.experiment.$i.fasta -r 0.10 -o seqs/hg38.experiment.$i.010.fasta
#done
python mutate_seqs.py -i seqs/hg38.random.fasta -r 0.01 -o seqs/hg38.random.001.fasta
python mutate_seqs.py -i seqs/hg38.random.fasta -r 0.02 -o seqs/hg38.random.002.fasta
python mutate_seqs.py -i seqs/hg38.random.fasta -r 0.03 -o seqs/hg38.random.003.fasta
python mutate_seqs.py -i seqs/hg38.random.fasta -r 0.04 -o seqs/hg38.random.004.fasta
python mutate_seqs.py -i seqs/hg38.random.fasta -r 0.05 -o seqs/hg38.random.005.fasta
python mutate_seqs.py -i seqs/hg38.random.fasta -r 0.06 -o seqs/hg38.random.006.fasta
python mutate_seqs.py -i seqs/hg38.random.fasta -r 0.07 -o seqs/hg38.random.007.fasta
python mutate_seqs.py -i seqs/hg38.random.fasta -r 0.08 -o seqs/hg38.random.008.fasta
python mutate_seqs.py -i seqs/hg38.random.fasta -r 0.09 -o seqs/hg38.random.009.fasta
python mutate_seqs.py -i seqs/hg38.random.fasta -r 0.10 -o seqs/hg38.random.010.fasta
