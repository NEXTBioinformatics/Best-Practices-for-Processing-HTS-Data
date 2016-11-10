Processing of DNA-seq
========

- Alignment
.. code-block:: bash
	/data/apps/bwa-0.7.12/bwa mem \
	-M \
	-t8 \
	-R"@RG\tID:exome_1\tPL:ILLUMINA\tPU:exome.11\tLB:$id\tSM:$id\tCN:AAU" \
	$assembly \
	$rawDataDir/$R1 \
	$rawDataDir/$R2 \
	> aligned_reads.sam

