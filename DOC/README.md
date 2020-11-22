Documentation of the functions: translate gffattr systime md5 edit_dist hamming end_adapter_pos charcount applytochars modstr setat

``systime()`` returns the number of milliseconds since the Linux epoch. This function is already in most other awk versions. Useful for timing.

``md5(str)`` returns the md5 code of the string argument. For example:
```
$ echo "example string to check" | bioawk_cas '{print; print "md5:", md5($0)}'
example string to check
md5: 59471d22e23e4198fd170cc7e4a58cbb
```
