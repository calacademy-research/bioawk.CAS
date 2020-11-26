Documentation of: translate gffattr systime md5 edit_dist hamming end_adapter_pos charcount applytochars modstr setat

``translate(nucl_str [, table_num])`` translates nucleotide string nucl_str and returns the protein sequence (from github.com/ctSkennerton/bioawk doc)

    bioawk -c fastx '{print ">"$name;print translate($seq)}' seq.fa.gz

can also use different translation tables. To translate using the bactera/archaea code:

    bioawk -c fastx '{print ">"$name;print translate($seq, 11)}' seq.fa.gz

``gffattr( attr_str, arr )`` parses the attr_str and puts the values in the arr argument. The attr_str is expected to be in the format of the gff attribute field.
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

``systime()`` returns the number of milliseconds since the Linux epoch. This function is already in most other awk versions. Useful for timing.

``md5(str)`` returns the md5 code of the string argument. For example:
```
$ echo "example string to check" | bioawk_cas '{print; print "md5:", md5($0)}'
example string to check
md5: 59471d22e23e4198fd170cc7e4a58cbb
```
A few of the functions were added since it is difficult or slow using substr() to modify the characters of a string.

``modstr`` takes 3 to 5 arguments ``modstr(str, start, length, [mod_type, str_length])`` and is used for in-place string variable case changing. mod_type 0 to lowercase, 1 to uppercase (default 0). Optional str_length faster for multiple calls, so the length isn't recalculated every call.

``setat`` takes 3 or 4 arguments ``setat(str,pos,replacement[, optional repeat_count])`` and does an in-place overwrite of a string with another string. Often used with a single character replacement string. The modified string is not changed in length.

``charcount(str, arr)`` fills arr with count of each character in str. returns number of different chars.

``applytochars(str, stmt_or_func)`` function calls the 2nd argument for each character in str with CHAR and ORD variables set.

``hamming( pattern, text [, text_pos: (1_indexed)default 1 [, case_sensitive: true [, N_wildcard: false] ]] )``
