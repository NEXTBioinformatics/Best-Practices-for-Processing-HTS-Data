Sketchy Pipeline

1. fastq files are imported to clc.
	Reads that did not pass a quality filter are ignored.
2. Trimming sequence ends.
	2a. Quality trimming: 
		Phred score > 20
		trim ambiguous nucleotides
	2b. Sequence filtering:
		Remove 1 3' terminal nucleotide
		number of residues of reads must be in that range: [30, 500]; short and long reads are discarded.
3. Map reads to reference (hg19). For algorithm see: http://www.clcbio.com/files/whitepapers/whitepaper-on-CLC-read-mapper.pdf
	Mapping parameters: match score 1, mismatch cost 2, linear gap cost, insertion cost 3, deletion cost 3, length fraction 0.95, similarity fraction 0.95, global alignment, ignore non-specific matches.
4. Local Realignment.
	Realign unaligned ends, perform local realignment 2 times.
5. Variant Detection. See: http://www.clcbio.com/files/whitepapers/whitepaper-probabilistic-variant-caller-1.pdf
	Minimum coverage:10, min count: 3, min frequency: 25.
	Restrict calling to target region: as given by SureSelect.