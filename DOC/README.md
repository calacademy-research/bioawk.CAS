#### Functions added to bioawk ####

``translate`` ``gffattr``  originally added in github.com/ctSkennerton/bioawk and ``bawk`` script added this repo is ``bioawk_cas -c fastx "$@"``

(1) ``translate(nucl_str [, table_num])`` translates nucleotide string nucl_str. Returns the protein sequence.

    bawk '{print ">"$name;print translate($seq)}' seq.fa.gz

can also use different translation tables. Here are the [genetic code table numbers](genetic_codes.md). To translate using the bacteria/archaea code table 11:

    bawk '{print ">"$name;print translate($seq, 11)}' seq.fa.gz

(2) ``gffattr( attr_str, arr )`` parses the attr_str and puts the values in the arr argument. The attr_str is expected to be in the format of the gff attribute field.
This field has subfields delimited by semi-colons where each subfield has a name and a value after an equal sign.  Returns number of subfields, which is same as length(arr).  For example, if this is the first line of the gff file:

```
NC_000085.6     Gnomon  gene    3064929 3075825 .       -       .       ID=gene-Gm38527;Dbxref=GeneID:102640308,MGI:MGI:5621412;Name=Gm38527;gbkey=Gene;gene=Gm38527;gene_biotype=lncRNA
```
```
bioawk 'NR==1{
   gffattr($NF, arr)
   print arr["ID"] " is the ID\n"
   for(subfield in arr)
      print subfield " = " arr[subfield]
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

```

Remaining functions have been added in bioawk_cas

**Miscellaneous functions** ``systime`` ``md5``

(3) ``systime()`` returns the number of milliseconds since the Linux epoch. This is already in most other awk versions. Useful for timing.

(4) ``md5(str)`` returns the md5 code of the string argument. For example:
```
$ echo "example string to check" | bioawk_cas '{print; print "md5:", md5($0)}'
example string to check
md5: 59471d22e23e4198fd170cc7e4a58cbb
```
**Character functions** ``modstr`` ``setat`` ``charcount`` ``applytochars``

A few of the functions were added since it is difficult or slow using substr() to modify the characters of a string.

(5) ``modstr`` takes 3 to 5 arguments ``modstr(str, start, length, [mod_type, str_length])`` and is used for in-place string variable case modification. mod_type 0 to lowercase, 1 to uppercase (default 0). Optional str_length faster for multiple calls, so the length isn't recalculated every call.

(6) ``setat`` takes 3 or 4 arguments ``setat(str,pos,replacement[, optional repeat_count])`` and does an in-place overwrite of a string with another string. Often used with a single character replacement string. The modified string is not changed in length.

(7) ``charcount(str, arr)`` fills arr with count of each character in str. returns number of different chars.

(8) ``applytochars(str, stmt_or_func)`` function calls the 2nd argument for each character in str with CHAR and ORD variables set.

**Search functions** ``hamming`` ``edit_dist`` ``end_adapter_pos``

(9) ``hamming( pattern, text [, text_pos: (1_indexed)default 1 [, case_sensitive: true [, N_wildcard: false] ]] )`` compares the pattern of the first argument to the characters in the text of the second argument for the length of the pattern (up to any remaining characters in the text string).

Returns the number of mismatches.

Comparison starts at position of the optional third argument, text_pos.  Default is to start at text string beginning which is text_pos 1.
Optional arguments 4 and 5 allow for case insensitive comparisons and the ability to treat the N character as a wildcard. To use either of these a text_pos must be provided.

(10) ``edit_dist( max_editdist, str1, str1_match_len, str2[, str2_len [, mode: default 1 [, flags]]] )``

           max_editdist: -1 means no max set. setting a max_editdist speeds up the search.   

           mode: 0 complete match, 1 prefix match, 2 infix match (add 10 or 20 for CIGAR). Can use string len -1 for full length.
           
           flags: 1 N matches ACTG, 2 Y matches CT, R matches AG, 3 both.
           
           Returns string starting with edit distance, -1 for no match, 0 for perfect and other values.


(11)  ``end_adapter_pos`` checks the last 16 nt of the sequence against the first 16 nt of the adapter seq, then the last 15 nt of the read for the first 15 nt of the adapter
and so-on until the last 4 nt of the read is checked with first 4 adapter nt. 

Considered matched when an attempt has an acceptable hamming distance: 4 mismatches at 16 nt, 3 starting at 12 nt then 2 starting at 8 nt, 1 mismatch at 6 and 5 nt, no mismatch at 4 nt.

``end_adapter_pos("", adapter)`` to set adapter, subsequently ``end_adapter_pos(seq)`` to check seq suffix against adapter prefix.

Returns string with 3 numbers: position of match, len, mismatches (-1 for none)
```
To set adapter call with empty seq and adapter as 2nd arg. subsequent calls use seq as only argument.
Only initial 16 bases of adapter are used.
       
end_adapter_pos("", "GATCGGAAGAGCACAC") # set to check first 16 bases of the TruSeq Adapter
end_adapter_pos($seq)                   # will check the last 16 bases of $seq against the set adapter
```    
For example
```
bawk 'BEGIN{end_adapter_pos("", "GATCGGAAGAGCACAC")}

   (rslt=end_adapter_pos($seq)) > 0 {
      print $name, "rec " NR":  ", rslt
      
 }' reads_L1_R1.fq.gz
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
 
 I have used this with libraries, especially mate-pair, after trimming with Trimmomatic to do an additional clean-up at the end of reads when the matched adapter is below the horizon of the other tool.
