MuTect2 Pitfalls
================

The standard parameters parameter in MuTect2 (GATK v3.7) are very strict. For example, a somatic variant (no matter frequency in the tumor sample) is discard if the variant is seen in more than one read in the normal sample, which is particularly problematic panels with read depth >1000X in targeted genes.

==================================== ======= ==========
Parameter                            Default Proposed
==================================== ======= ==========
max_alt_alleles_in_normal_count      1       10000000
max_alt_allele_in_normal_fraction    0.03    0.10
==================================== ======= ==========

We propose solving this issue by raising the values and adding downstream filters (discussed in sections below) to filter out the inevitable false positives. By setting `max_alt_alleles_in_normal_count` to a very high number, we never filter out variants based on the absolute count of ALT alleles in the normal sample, and we allow up to 10% of the normal reads to contain the ALT allele.

Motivation
----------

Coming soon...

Downstream Filters
------------------

Coming soon...
