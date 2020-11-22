Documentation of the functions: translate gffattr systime md5 edit_dist hamming end_adapter_pos charcount applytochars modstr setat

``systime()`` returns the number of milliseconds since the Linux epoch. This function is already in most other awk versions. Useful for timing.

``md5(str)`` returns the md5 code of the string argument. For example:
```
$ echo "example string to check" | bioawk_cas '{print; print "md5:", md5($0)}'
example string to check
md5: 59471d22e23e4198fd170cc7e4a58cbb
```
A few of the functions were added since it is difficult or slow using substr() to modify the characters of a string.

``modstr`` takes 3 to 5 arguments ``modstr(str, start, length, [mod_type, str_length])`` and is used for in-place string variable case changing. mod_type 0 to lowercase, 1 to uppercase (default 0). Optional str_length faster for multiple calls, so the length isn't recalculated every call.

``setat`` takes 3 or 4 arguments ``setat(str,pos,replacement[, optional repeat_count])`` and does an in-place overwrite of a string with another string. Often used with a single character replacement string. The modified string in not changed in length.

``charcount(str, arr)`` fills arr with count of each character in str. returns number of different chars.

``applytochars(str, stmt_or_func)`` function calls the 2nd argument called for each character in str with CHAR and ORD variables set.
