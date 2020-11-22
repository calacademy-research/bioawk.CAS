#ifndef ADDON_H
#define ADDON_H

/* do not change these values, as they are tied to col_def[][] in addon.c */
#define BIO_NULL -1
#define BIO_HDR   0
#define BIO_BED   1
#define BIO_SAM   2
#define BIO_VCF   3
#define BIO_GFF   4
#define BIO_FASTX 5

#define BIO_SHOW_HDR 0x1

extern int bio_fmt, bio_flag;
extern char *bio_hdr_chr;

int bio_get_fmt(const char *s);
int bio_skip_hdr(const char *r);
void bio_set_colnm(void);

int bio_getrec(char **pbuf, int *psize, int isrecord);

/* The following explains how to add a new function. 1) Add a function index
 * (e.g. #define BIO_FFOO 102) in addon.h. The integer index must be larger than
 * 14 in the current awk implementation (see also macros starting with "F"
 * defined in awk.h). 2) Add the function name and the function index in the
 * keywords array defined in lex.c. Remember to keep the array sorted. 3)
 * Implement the actual function in bio_func(). */

#define BIO_FAND      101
#define BIO_FOR       102
#define BIO_FXOR      103
#define BIO_FLSHIFT   104 /* to add */
#define BIO_FRSHIFT   105 /* to add */
#define BIO_FCOMPL    106 /* to add */

#define BIO_FREVERSE  201
#define BIO_FREVCOMP  202
#define BIO_FGC       203
#define BIO_FMEANQUAL 204
#define BIO_FQUALCOUNT 205
#define BIO_FTRIMQ    206
#define BIO_TRANSLATE 207
#define BIO_FSETAT    208 /* directly set chars at a position (index+1) of string variable JBH 06Mar2019 */
#define BIO_FSYSTIME  209 /* implement systime() which is an gawk extension JBH 12Aug2019 */
#define BIO_FHAMMING  210 /* hamming(pat, text, compare_len[, optional case_sensitive]) */
#define BIO_FEDLIB    211 /* entry point for calling edlib JBH 29Mar2020 */
#define BIO_ADAPATEND 212 /* end_adapter_pos looks for an adapter prefix at sequence end */
#define BIO_FMD5      213 /* return MD5 string of arg1 in md5(str_arg) */
#define BIO_CHARCOUNT 214 /* charcount(str, ar_chars) -- returns val e.g. ar_chars["A"]=173 JBH 21Jul2020 */
#define BIO_MODSTR    215 /* modstr(str, start, length, mod_type) -- mod_type 0 to lowercase, 1 to uppercase, might add ability to use function for mod_type later */
#define BIO_APPLYCHARS 216 /* applytochars(str, eval_for_each_char_in_str[, eval_for_each_char...]) -- loops over chars in str executing arg with var CHAR set */
#define BIO_GFFATTR   217 /* added in github.com/ctSkennerton/bioawk */

struct Cell;
struct Node;

struct Cell *bio_func(int f, struct Cell *x, struct Node **a);

#endif
