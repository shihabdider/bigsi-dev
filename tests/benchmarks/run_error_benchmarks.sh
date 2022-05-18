for i in {5..20};
do
	time python benchmark.py -q seqs/hg38.experiment.$i.001.fasta -c hg38.office.config.json -o output/hg38.experiment.$i.001 -i 99
	time python benchmark.py -q seqs/hg38.experiment.$i.002.fasta -c hg38.office.config.json -o output/hg38.experiment.$i.002 -i 98
	time python benchmark.py -q seqs/hg38.experiment.$i.003.fasta -c hg38.office.config.json -o output/hg38.experiment.$i.003 -i 97
	time python benchmark.py -q seqs/hg38.experiment.$i.004.fasta -c hg38.office.config.json -o output/hg38.experiment.$i.004 -i 96
	time python benchmark.py -q seqs/hg38.experiment.$i.005.fasta -c hg38.office.config.json -o output/hg38.experiment.$i.005 -i 95
	time python benchmark.py -q seqs/hg38.experiment.$i.006.fasta -c hg38.office.config.json -o output/hg38.experiment.$i.006 -i 94
	time python benchmark.py -q seqs/hg38.experiment.$i.007.fasta -c hg38.office.config.json -o output/hg38.experiment.$i.007 -i 93
	time python benchmark.py -q seqs/hg38.experiment.$i.008.fasta -c hg38.office.config.json -o output/hg38.experiment.$i.008 -i 92
	time python benchmark.py -q seqs/hg38.experiment.$i.009.fasta -c hg38.office.config.json -o output/hg38.experiment.$i.009 -i 91
	time python benchmark.py -q seqs/hg38.experiment.$i.010.fasta -c hg38.office.config.json -o output/hg38.experiment.$i.010 -i 90
done
#time python benchmark.py -q seqs/hg38.random.007.fasta -c hg38.office.config.json -o output/hg38.random.007 -i 93
#time python benchmark.py -q seqs/hg38.random.008.fasta -c hg38.office.config.json -o output/hg38.random.008 -i 92
#time python benchmark.py -q seqs/hg38.random.009.fasta -c hg38.office.config.json -o output/hg38.random.009 -i 91
