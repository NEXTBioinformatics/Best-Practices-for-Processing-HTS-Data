Software and reference data
========
- Quality check of raw reads
	Quality check and removal of adapters from raw reads is done using the wrapper tool "Trim Galore!" which combines adapter removal with Cutadapt and quality checks with FastQC:
		- `Trim Galore! <https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/>`_
		- `Cutadapt <https://cutadapt.readthedocs.io/en/stable/>`_
		- `FastQC <https://www.bioinformatics.babraham.ac.uk/projects/fastqc/>`_

- GATK best practices
	These workflows are based on the `GATK Best Practices <https://software.broadinstitute.org/gatk/best-practices/>`_, with the addition of a second variant caller. The workflow requires specification of paths to a
	number of programs and reference datasets which must be downloaded and installed first:
		- `BWA <http://bio-bwa.sourceforge.net/>`_
		- `GATK <https://software.broadinstitute.org/gatk/download/>`_
		- `Picard Tools <http://broadinstitute.github.io/picard/>`_
		- `SAMtools <http://www.htslib.org/>`_
		- `VarScan2 <http://varscan.sourceforge.net/index.html>`_

- Reference genome
	In the NEXT bioinformatics network we use the Genomic Data Commons (GDC) version of the GRCh38 reference genome. 
	the reference genome and associated index files is available for download `here
	<https://gdc.cancer.gov/about-data/data-harmonization-and-generation/gdc-reference-files>`_
	
- GATK bundle
	For the GATK workflow a number of reference datasets with known variants are needed. A `ressource bundle
	<ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/>`_ 
	with all necessary files for the GATK workflow is provided by the Broad Institute.
	
- COSMIC
	Somatic variant calling using Mutect2 uses a whitelist of mutations previously seen in cancers saved in a VCF file. 
	VCF files containing both coding and non-coding mutations can be downloaded by following instructions on the COSMIC `download page <http://cancer.sanger.ac.uk/cosmic/download>`_.