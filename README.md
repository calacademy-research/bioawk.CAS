**bioawk_cas** is a fork of Connor Skennerton's fork https://github.com/ctSkennerton/bioawk of Heng Li's https://github.com/lh3/bioawk. Connor adds a translate function to the core bioawk extensions to be able to translate nucelotide into protein sequences.

This version, bioawk_cas, adds several additional functions including linking in edlib to perform approximate searches. The original documention is below these notes.

To see the novel bioawk functions use -h or --help,
```
bioawk_cas -h

usage: bioawk_cas [-F fs] [-v var=value] [-c fmt] [-tH] [-f progfile | 'prog'] [file ...]

bed:
	1:chrom 2:start 3:end 4:name 5:score 6:strand 7:thickstart 8:thickend 9:rgb 10:blockcount 11:blocksizes 12:blockstarts 
sam:
	1:qname 2:flag 3:rname 4:pos 5:mapq 6:cigar 7:rnext 8:pnext 9:tlen 10:seq 11:qual 
vcf:
	1:chrom 2:pos 3:id 4:ref 5:alt 6:qual 7:filter 8:info 
gff:
	1:seqname 2:source 3:feature 4:start 5:end 6:score 7:filter 8:strand 9:group 10:attribute 
fastx:
	1:name 2:seq 3:qual 4:comment 

bioawk functions:
	gc meanqual qualcount revcomp reverse trimq and or xor
	translate systime md5 edit_dist hamming end_adapter_pos charcount applytochars modstr setat
```
The first line under bioawk functions in the above code block are the functions added in Heng Li's original version.
The next line has the translate function and then the, currently 9, new functions added in bioawk_cas.

Since the most common use of bioawk is with fasta or fastq files using the -c fastx option, a script named **bawk** is included that presumes this.
**bawk** is `bioawk -c fastx "$@"` and saves a bit of typing.

Most of the new functions will provide a parameter overview if you use its name in a BEGIN block. For example,
```
bioawk_cas 'BEGIN{md5}'
bioawk_cas: md5 takes a string argument and returns its md5 value.
```
or
```
bawk 'BEGIN{edit_dist}'
bioawk_cas: edit_dist requires 4 to 7 arguments: max_editdist, str1, str1_match_len, str2[, str2_len [, mode: default 1 [, flags]]]
                  mode: 0 complete match, 1 prefix match, 2 infix match (add 10 or 20 for CIGAR). Can use string len -1 for full length.
                  flags: 1 N matches ACTG, 2 Y matches CT, R matches AG, 3 both.
```
An edlib object file for Linux 86_64 systems and one for macOS can be used by the make file if edlib.obj does not exist.
If these do not work, clone the edlib repo https://github.com/Martinsos/edlib and after running make copy the edlib.obj into the bioawk_cas repo directory.

Examples and function documentation to come.

**Here is the original documentation**

###Introduction

Bioawk is an extension to [Brian Kernighan's awk][1], adding the support of
several common biological data formats, including optionally gzip'ed BED, GFF,
SAM, VCF, FASTA/Q and TAB-delimited formats with column names. It also adds a
few built-in functions and an command line option to use TAB as the
input/output delimiter. When the new functionality is not used, bioawk is
intended to behave exactly the same as the original BWK awk.

###New functionality

#####Command line option `-t` *note: tab character for FS and OFS is now the default*

Using this option is equivalent to

    bioawk -F'\t' -v OFS="\t"

#####Command line option `-c arg`

This option specifies the input format. When this option is in use, bioawk will
seamlessly add variables that name the fields, based on either the format or
the first line of the input, depending *arg*. This option also enables bioawk
to read gzip'd files. The argument *arg* may take the following values:

* `help`. List the supported formats and the naming variables.

* `hdr` or `header`. Name each column based on the first line in the input.
  Special characters in the first are converted to underscore. For example:

        grep -v ^## in.vcf | bioawk -tc hdr '{print $_CHROM,$POS}'

  prints the `CHROM` and `POS` columns of the input VCF file.

* `sam`, `vcf`, `bed` and `gff`. SAM, VCF, BED and GFF formats.

* `fastx`. This option regards a FASTA or FASTQ as a TAB delimited file with
  four columns: sequence name, sequence, quality and FASTA/Q comment, such that
  various fields can be retrieved with column names. See also example 4 in the
  following.

#####New built-in functions

See `awk.1`.

###Examples

1. List the supported formats:

        bioawk -c help

2. Extract unmapped reads without header:

        bioawk -c sam 'and($flag,4)' aln.sam.gz

3. Extract mapped reads with header:

        bioawk -Hc sam '!and($flag,4)'

4. Reverse complement FASTA:

        bioawk -c fastx '{print ">"$name;print revcomp($seq)}' seq.fa.gz

5. Create FASTA from SAM (uses revcomp if FLAG & 16)

        samtools view aln.bam | \
            bioawk -c sam '{s=$seq; if(and($flag, 16)) {s=revcomp($seq)} print ">"$qname"\n"s}'

6. Print the genotypes of sample `foo` and `bar` from a VCF:

        grep -v ^## in.vcf | bioawk -tc hdr '{print $foo,$bar}'

7. Translate nucleotide into protein sequence
 
        bioawk -c fastx '{print ">"$name;print translate($seq)}' seq.fa.gz
can also use different translation tables.  To translate using the
bactera/archaea code:

        bioawk -c fastx '{print ">"$name;print translate($seq, 11)}' seq.fa.gz



###Potential limitations

1. When option `-c` is in use, bioawk replaces the line reading module of awk.
   The new line reading function parses FASTA and FASTQ files and seamlessly
   reads gzip'ed files. However, the new code does not fully mimic the original
   code. It may fail in corner cases (though this has not happened yet). Thus
   when `-c` is not specified, awk falls back to the original line reading code
   and does not support gzip'ed input.

2. When `-c` is in use, several strings allocated in the new line reading
   module are not freed in the end. These will be reported by valgrind as
   "still reachable". To some extent, these are not memory leaks.


[1]: http://www.cs.princeton.edu/~bwk/btl.mirror/
