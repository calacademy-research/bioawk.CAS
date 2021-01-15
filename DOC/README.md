### Functions added to bioawk ###

**ctSkennerton's and derivative functions** ``translate`` ``gffattr`` ``gtfattr`` ``samattr``

``translate`` ``gffattr``  originally added in https://github.com/ctSkennerton/bioawk and ``bawk`` script added here is ``bioawk_cas -c fastx "$@"``. Optional arg added to this version of gffattr. Also ``gtfattr`` to parse the slightly different syntax of gtf attributes and ``samattr`` for sam tags.

(1) ``translate(nucl_str [, table_num])`` translates nucleotide string nucl_str. Returns the protein sequence.

    bawk '{print ">"$name;print translate($seq)}' seq.fa.gz

can also use different translation tables. Here are the [genetic code table numbers](genetic_codes.md). To translate using the bacteria/archaea code table 11:

    bawk '{print ">"$name;print translate($seq, 11)}' seq.fa.gz

(2) ``gffattr( attr_str, arr[, pos_arr])`` or ``gtfattr( attr_str, arr[, pos_arr])`` parses the attr_str and puts the values in the arr argument. The attr_str is expected to be in the format of the gff attribute field (gtf attribute field format for gtfattr).
This field has subfields delimited by semi-colons where each subfield has a name and a value after an equal sign (space for gtfattr). Quotes removed around value for gtfattr.

Returns number of subfields, which is same as length(arr). Optional pos_array will put the first key name in pos_arr[1], second in pos_arr[2], etc.
The pos_arr option added in bioawk_cas.
For example, if this is the first line of the gff file:

```
NC_000085.6     Gnomon  gene    3064929 3075825 .       -       .       ID=gene-Gm38527;Dbxref=GeneID:102640308,MGI:MGI:5621412;Name=Gm38527;gbkey=Gene;gene=Gm38527;gene_biotype=lncRNA
```
```
bioawk_cas 'NR==1{
   gffattr($NF, arr)
   print arr["ID"] " is the ID\n"
   
   for(subfield in arr)
      print subfield " = " arr[subfield]
      
   print "\nOutput in same order as attr string"
   flds = gffattr($NF, arr, pos)
   for(f = 1; f <= flds; f++)
      printf "%d\t%s\t%s\n", f, pos[f], arr[ pos[f] ]
}' example.gff
```
outputs
```
gene-Gm38527 is the ID

gbkey = Gene
ID = gene-Gm38527
gene_biotype = lncRNA
Dbxref = GeneID:102640308,MGI:MGI:5621412
gene = Gm38527
Name = Gm38527

Output in same order as attr string
1	ID	gene-Gm38527
2	Dbxref	GeneID:102640308,MGI:MGI:5621412
3	Name	Gm38527
4	gbkey	Gene
5	gene	Gm38527
6	gene_biotype	lncRNA
```

Remaining functions have been added in bioawk_cas

(3) ``samattr(sam_line, arr[, pos_arr])`` similar to gffattr and gtfattr except since tag fields are tab delimited, the entire line is provided to the function using the $0 variable. The type char is also in pos_arr with 'T' prefixed position, pos_arr["T1"], pos_arr["T2"], etc. For example, if this was the sam line:
```
ref1_grp1_p004  99      ref1    13      6       10M     =       37      34      CCGGGGATCC      ''''''''''      fa:f:1.38e-23   za:Z:xRG:Z:grp2 RG:Z:grp1       NM:i:0  MD:Z:10
```
then
```
bioawk_cas 'BEGIN{OFS="\t"}
!/^@/ {
    tot = samattr($0, ar, pos)
    for(t=1; t<=tot; t++) {
      key = pos[t]; val = ar[key]
      typ = pos["T" t]
      print t, key,  val, typ
    }
} ' oneliner_example.sam
```
outputs
```
1       fa      1.38e-23    f
2       za      xRG:Z:grp2  Z
3       RG      grp1    Z
4       NM      0   i
5       MD      10  Z
```

**Miscellaneous functions** ``systime`` ``md5``

(4) ``systime()`` returns the number of milliseconds since the Linux epoch. This is already in most other awk versions. Useful for timing.

(5) ``md5(str)`` returns the md5 code of the string argument. For example:
```
$ echo "example string to check" | bioawk_cas '{print; print "md5:", md5($0)}'
example string to check
md5: 59471d22e23e4198fd170cc7e4a58cbb
```
(6) ``fldcat(start_fldno, end_fldno[, separator])`` return fields from start_fldno to end_fldno with separator given or OFS if not. This is useful to print the fields starting from one field to the end of fields using NF for end_fldno. Or, for example if you could give just 12th to final field to samattr instead of $0.
```
flds = fldcat(12, NF, "\t")
tot=samattr(flds, ar, pos)
...
```

**Character functions** ``modstr`` ``setat`` ``charcount`` ``applytochars``

A few of the functions were added since it is difficult or slow using substr() to modify or access the characters of a string.

(7) ``modstr`` takes 3 to 5 arguments ``modstr(str, start, length, [mod_type, str_length])`` and is used for in-place string variable case modification. mod_type 0 to lowercase, 1 to uppercase (default 0). Optional str_length faster for multiple calls, so the length isn't recalculated every call.

(8) ``setat`` takes 3 or 4 arguments ``setat(str,pos,replacement[, optional repeat_count])`` and does an in-place overwrite of a string with another string. Often used with a single character replacement string. The modified string is not changed in length.

(9) ``charcount(str, arr)`` fills arr with count of each character in str. returns number of different chars.
```
bioawk_cas 'BEGIN{OFS=":"
   exmp="AaGBCNdEfaGHNINNJ"
   charcount(exmp, arr)
   for(c in arr)
      print c, arr[c]
}' |  sort -t ":" -k2,2nr | tr "\n" " "

```
gives
```
N:4 G:2 a:2 A:1 B:1 C:1 E:1 H:1 I:1 J:1 d:1 f:1
```

(10) ``applytochars(str, stmt_or_func [,...])`` call the 2nd and other arguments for each character in str with CHAR and ORD variables set.
```
bioawk_cas '
   function apply(){if(CHAR=="N")Ns++}
   function prt(s){print s}
   BEGIN{
      str="AaNBCNdEfNGHNINNJ"

      # applytochars(str, apply())  # example just using the function to count the Ns
      applytochars(str, Ns += CHAR=="N", prt(CHAR ":" ORD ":" or(ORD, 32))) # example using a statement and a function

      # CHAR and ORD still set to last char in string
      print "\nLast CHAR is:", CHAR, ORD, "or32", or(ORD, 32)
      print "tot Ns:", Ns
}'
```
gives
```
A:65:97
a:97:97
N:78:110
B:66:98
C:67:99
N:78:110
d:100:100
E:69:101
f:102:102
N:78:110
G:71:103
H:72:104
N:78:110
I:73:105
N:78:110
N:78:110
J:74:106

Last CHAR is: J 74 or32 106
tot Ns: 6
```

**Search functions** ``hamming`` ``edit_dist`` ``end_adapter_pos``

(11) ``hamming( pattern, text [, text_pos: (1_indexed)default 1 [, case_sensitive: true [, N_wildcard: false] ]] )`` compares the pattern of the first argument to the characters in the text of the second argument for the length of the pattern (up to any remaining characters in the text string).

Returns the number of mismatches.

Comparison starts at position of the optional third argument, text_pos.  Default is to start at text string beginning which is text_pos 1.
Optional arguments 4 and 5 allow for case insensitive comparisons and the ability to treat the N character as a wildcard. To use either of these a text_pos must be provided.

(12) ``edit_dist( max_editdist, str1, str1_match_len, str2[, str2_len [, mode: default 1 [, flags]]] )``

           max_editdist: -1 means no max set. setting a max_editdist speeds up the search.   

           mode: 0 complete match, 1 prefix match, 2 infix match (add 10 or 20 for CIGAR). Can use string len -1 for full length.
           
           flags: 1 N matches ACTG; 2 Y matches CT, R matches AG; 3 both.
           
           Returns string starting with edit distance, -1 for no match 0 for exact match, and other values.

Example below looks for a PacBio amplification adapter in CCS reads.
Search mode is 22: 2 for InFix (HW) search plus 20 to add the extended CIGAR match to the output.
The edit_dist() call is in the adapterMatch function which is called for forward and for reverse compliment of the adapter.
```
bawk '
   function int_ceil(f) { add = !(f==int(f)); return int(f)+add }
   function adapterMatch(Adapter, AdapStr) {
      matchInfo = edit_dist(max_miss, Adapter, alen, $seq, slen, mode)
      if(matchInfo >= 0) {
         printf "%s %s %s", prefix, AdapStr, matchInfo
         prefix = "\t"
      }
   }

   BEGIN{InFix = 2; ExtCIGAR = 20; RegCIGAR = 10; extmode = InFix + ExtCIGAR; ComputeLen = -1

         adap="AAGCAGTGGTATCAACGCAGAGTACT"; adap_rc=adap; revcomp(adap_rc)
         alen=length(adap)
         max_miss = int_ceil(alen/10) + 1 # little less than 90% match at worst

         #mode = InFix
         mode = extmode
   }
   # check for amplification adapter
   {  slen = length($seq)
      prefix = $name "\t"

      adapterMatch(adap,    "ampF")
      adapterMatch(adap_rc, "ampR")

      if(prefix == "\t")
         print "\t rdlen " slen
 }' subreads_ccs.fasta
```
This outputs in part
```
m64044_201011_075919/1/ccs       ampF 0 26= 1 26         ampR 0 26= 11994 12019  rdlen 12019
m64044_201011_075919/2/ccs       ampF 0 26= 1 26         ampR 0 26= 17194 17219  rdlen 17219
m64044_201011_075919/3/ccs       ampF 0 26= 1 26         ampR 2 12=1D7=1D7= 14871 14898  rdlen 14898
m64044_201011_075919/7/ccs       ampF 1 5=1I20= 1 25     ampR 3 8=1I4=1X5=1X6= 19327 19351       rdlen 19351
m64044_201011_075919/8/ccs       ampF 0 26= 1 26         ampR 1 25=1I 10536 10560        rdlen 10560
m64044_201011_075919/48/ccs      ampF 0 26= 1 26         ampR 0 26= 11256 11281  rdlen 11281
m64044_201011_075919/49/ccs      ampF 1 5=1I20= 1 25     ampR 2 19=1X5=1I 8912 8936 8912 8937 8912 8938  rdlen 8938
m64044_201011_075919/51/ccs      ampF 1 10=1D16= 1 27    ampR 0 26= 11678 11703  rdlen 11703
m64044_201011_075919/118/ccs     ampF 2 17=1X5=1D3= 1 27         ampR 0 26= 12642 12667  rdlen 12667
m64044_201011_075919/120/ccs     ampF 0 26= 1 26         ampR 1 4=1D22= 17923 17949      rdlen 17949
m64044_201011_075919/124/ccs     ampF 0 26= 1 26         ampR 0 26= 13041 13066  rdlen 13066
```

(13)  ``end_adapter_pos`` checks the last 16 nt of the sequence against the first 16 nt of the adapter seq, then the last 15 nt of the read for the first 15 nt of the adapter
and so-on until the last 4 nt of the read is checked with first 4 adapter nt. 

Considered matched when an attempt has an acceptable hamming distance: 4 mismatches at 16 nt, 3 starting at 12 nt then 2 starting at 8 nt, 1 mismatch at 6 and 5 nt, no mismatch at 4 nt.

``end_adapter_pos("", adapter)`` to set adapter, subsequently ``end_adapter_pos(seq)`` to check seq suffix against adapter prefix.

Returns string with 3 numbers: position of match, len, mismatches (-1 for none)
    
For example
```
bawk 'BEGIN{end_adapter_pos("", "GATCGGAAGAGCACAC")} # set to check first 16 bases of the TruSeq Indexed Adapter

   (rslt=end_adapter_pos($seq)) > 0 {  # check the last 16 to last 4 bases of $seq against the adapter head for close match
   
      split(rslt, ar, " ")
      pos = ar[1]; mlen = ar[2]; mismatches = ar[3]
      
      printf "%s\trec %s:  \t%s %s %s\n", $name, NR, pos, mlen, mismatches
      
 }' reads_L4_R1.fq.gz
 ```
 outputs in part
 ```
E00489:558:H7N5LCCX2:4:1101:24718:1344  rec 94:         139 11 3
E00489:558:H7N5LCCX2:4:1101:27600:1344  rec 103:        143 7 0
E00489:558:H7N5LCCX2:4:1101:23531:1397  rec 165:        146 4 0
E00489:558:H7N5LCCX2:4:1101:26951:1415  rec 189:        134 16 3
E00489:558:H7N5LCCX2:4:1101:23429:1432  rec 217:        135 15 0
E00489:558:H7N5LCCX2:4:1101:25733:1450  rec 245:        143 7 2
E00489:558:H7N5LCCX2:4:1101:24698:1485  rec 287:        145 5 0
 ```
 Top line tells us record 94 has the first 11 bases of the adapter matching with 3 errors starting at position 139 of the read.
 
 I have used this with libraries, especially mate-pair, after trimming with Trimmomatic to do an additional clean-up at the end of reads when the matched adapter length is below the horizon of the other tool.
