Find potential de novo variants in a given VCF.

de_novo_finder_3.py is a major update to the original caller. The main purpose is to
identify events that appear to be de novo in a specified VCF file that contains
sequence information from trios. Due to the nature of this variation, parents are required
to be homozygous reference, while children should usually be heterozygous for the
variant. Confidence in the calls is established with quality filters:

1) The variant must pass all of the filters applied by the variant caller
    To accept TruthSensitivityTranche variants, use the -q flag.
    ** Removed in 3.93 **
2) The PL (normalized, Phred-scaled likelihoods for AA, AB, and BB genotypes where
A = ref and B = alt) score of the child is required to be >T : 0 : >0 for a given
threshold, T.
    T is set at a default of 20. It can be adjusted with the -t flag.
3) The allelic balance (# alternative reads/total reads) of the child is required
to be at least 20%. The allelic balance in the parents should be less than or
equal to 5%.
    These numbers can be adjusted with the -c and -p flags, respectively.
4) The depth in the child is required to be greater than a tenth of the sum of the
depths in both parents.
    The fraction of the sum of depths can be adjusted with the -d flag.

This script processes both single nucleotide variants (SNVs) and small insertions and
deletions (indels). To skip indels, use the -i flag. ** Removed in 3.81 **

Lines in the VCF that have multiple alternative alleles are processed only if all
alleles are single bases.

Instead of requiring a hard PL threshold in the parents, we have defined a relative
probability of an event being truly de novo versus the probability that it was a missed
heterozygote call in one of the two parents (the most likely error mode).
p_dn = P(true de novo | data) / (P(true de novo | data) + P(missed het in parent | data))

where P(true de novo | data) = P(data | true de novo) * P(true de novo)
P(data | true de novo) = Pdad_ref * Pmom_ref * Pchild_het
P(true de novo) = 1/30 Mb
and P(missed het in parent | data) = P(data | at least one parent is het) * P(one parent het)
P(data | at least one parent is het) = (Pdad_ref*Pmom_het + Pdad_het*Pmom_ref) * Pchild_het
P(one parent het) = 1 - (1-f)^4
where f is the maximum of the frequency of the variant in the VCF or ESP

The minimum p_dn considered is 0.05, but can be adjusted with the -m flag.

The potential de novo variants are then split by SNVs and indels into HIGH, MEDIUM,
and LOW validation likelihood by the following criteria.
HIGH_SNV:
p_dn > 0.99 and child_AD > 0.3 and dp_ratio (kid depth/combined parental depth) > 0.2
or
p_dn > 0.99 and child_AD > 0.3 and allele count (AC) == 1

MEDIUM_SNV:
p_dn > 0.5 and child_AD > 0.3
or
p_dn > 0.5 and child_AD  > 0.2 and AC == 1

LOW_SNV:
p_dn > 0.05 and child_AD > 0.2


HIGH_indel:
p_dn > 0.99 and child_AD > 0.3 and dp_ratio > 0.2
or
p_dn > 0.99 and child_AD > 0.3 and AC == 1

MEDIUM_indel:
p_dn > 0.5 and child_AD > 0.3
or
p_dn > 0.5 and child_AD > 0.2 and AC ==1

LOW_indel:
p_dn > 0.05 and child_AD > 0.2


If SnpEff annotations have been included in the annotation line of the VCF, the -a
flag can be used to extract and print the gene name and mutation category. The same is
true for VEP annotations using the -v flag.

The updates to this caller require changes in the input. A PED file is required to
establish the family relations. PED format has 6 columns: family ID, child, dad, mom,
sex of the child, affected status of the child. The ESP counts file is required to
determine the population frequency of an event.

Current ESP counts file: all_ESP_counts_5.28.13.txt

Output contains many more columns than earlier versions of the script so downstream
scripts may need to be adjusted.


3.4: different flags
3.5: fixed depth bug and now prints DP_ratio, child sex and affected status
3.52: modified to remove "chr" if it is in the "chr" column and skip lines
that have no PL information
3.6: Changed the INCORRECT MEDIUM_indel variant AC requirement. It used
to be AC >= 5, but it should really be AC <= 5. (This may been seen in
some versions of 3.52.) In addtion, 3.6 skips lines with no AD information
and keeps track of the number.
3.7: Adjusted validation likelihood filters so that MEDIUM_SNV that meet
the following criteria are moved to HIGH_SNV:
    AC < 10
    child AD > 0.3
    child depth > 10
3.71: Added ability to handle .gz VCFs
3.72: Haplotype caller for 0/0 individuals will list an AD of "." which breaks the
script. I am now assuming that there are 0 alternative reads in these cases.
3.73: Haplotype caller now for the hemizygous variants
3.74: Fixed the labels reading line so that it strips off the newline character
3.75: Added the ability to handle VEP annotations (curtesy of Jack Kosmicki)
3.8: The VEP annotation section was completely wrong. Rewrote with code from Konrad
Karczewski loftee_utils.py (https://github.com/konradjk/loftee/blob/master/src/loftee_utils.py)
3.9: DROPPED THE -i FLAG. The script does not throw out multiallelic lines where a SNV
and indel are present. It, however, still does not handle multialleleic lines with more
than two alternative alleles.
3.91: Slight update to VEP annotation list
3.92: Adding a quick fix to avoid situations where SnpEff where EFFECT (such as
'INTRON') is missing
3.93: Removed the -q flag that allowed individuals to look at Truth Sensitivity
Tranche variants
3.94: Had never incorporated that if the greatest frequency (f) == 0, then it
becomes 100/30Mbp
3.95: Fixed a bug that led to the first variant being dropped from the VCF
