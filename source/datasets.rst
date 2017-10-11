Software and reference data
========
- GATK best practices
	These workflows are based on the `GATK Best Practices <https://software.broadinstitute.org/gatk/best-practices/>`_. The workflow requires specification of paths to a
	number of programs and reference datasets which must be downloaded and installed first:
	`BWA <http://bio-bwa.sourceforge.net/>`_
	`GATK <https://software.broadinstitute.org/gatk/download/>`_
	`Picard Tools <http://broadinstitute.github.io/picard/>`_
	`SAMtools <http://www.htslib.org/>`_
	`VarScan2 <http://varscan.sourceforge.net/index.html>`_

- Reference genome
	In the NEXT bioinformatics network we use the Genomic Data Commons (GDC) version of the GRCh38 reference genome. 
	The reference genome and associated index files is available for download `here
	<https://gdc.cancer.gov/about-data/data-harmonization-and-generation/gdc-reference-files>`_
	
- GATK bundle
	For the GATK workflow a number of reference datasets with known variants are needed. A `ressource bundle
	<ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/>`_ 
	with all necessary files for the GATK workflow is provided by the Broad Institute.