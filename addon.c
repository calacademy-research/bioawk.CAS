#include <math.h>
#include <ctype.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include "awk.h"
#include "edlib.h"
#include "end_adapter.h"
extern char *md5str(uint8_t *msg, size_t len);

int bio_flag = 0, bio_fmt = BIO_NULL;

static const char *col_defs[][15] = { /* FIXME: this is convenient, but not memory efficient. Shouldn't matter. */
    {"header", NULL},
    {"bed", "chrom", "start", "end", "name", "score", "strand", "thickstart", "thickend", "rgb", "blockcount", "blocksizes", "blockstarts", NULL},
    {"sam", "qname", "flag", "rname", "pos", "mapq", "cigar", "rnext", "pnext", "tlen", "seq", "qual", NULL},
    {"vcf", "chrom", "pos", "id", "ref", "alt", "qual", "filter", "info", NULL},
    {"gff", "seqname", "source", "feature", "start", "end", "score", "filter", "strand", "group", "attribute", NULL},
    {"fastx", "name", "seq", "qual", "comment", NULL},
    {NULL}
};

static const char *tab_delim = "nyyyyyn", *hdr_chr = "\0#@##\0\0";

/************************
 * Setting column names *
 ************************/

static void set_colnm_aux(const char *p, int col)
{
    const char *q;
    char *r = 0;
    Cell *x;
    for (q = p; *q; ++q) /* test if there are punctuations */
        if (ispunct(*q) && *q != '_') break;
    if (*q || isdigit(*p)) { /* there are punctuations or the first is digit */
        char *qq;
        r = malloc(strlen(p) + 2);
        if (isdigit(*p)) {
            *r = '_';
            strcpy(r + 1, p);
        } else strcpy(r, p);
        for (qq = r; *qq; ++qq)
            if (ispunct(*qq)) *qq = '_';
        q = r;
    } else q = p;
    if ((x = lookup(q, symtab)) != NULL) /* do not add if not appear in the program */
        setfval(x, (Awkfloat)col);
    if (r) free(r);
}

int bio_get_fmt(const char *s)
{
    int i, j;
    if (strcmp(s, "hdr") == 0) return BIO_HDR;
    for (i = 0; col_defs[i][0]; ++i)
        if (strcmp(s, col_defs[i][0]) == 0) return i;
    for (i = 1; col_defs[i][0]; ++i) {
        printf("%s:\n\t", col_defs[i][0]);
        for (j = 1; col_defs[i][j]; ++j)
            printf("%d:%s ", j, col_defs[i][j]);
        putchar('\n');
    }
    return BIO_NULL;
}

int bio_skip_hdr(const char *r)
{
    if (bio_fmt <= BIO_HDR) return 0;
    if (*r && *r == hdr_chr[bio_fmt]) {
        if (bio_flag & BIO_SHOW_HDR) puts(r);
        return 1;
    } else return 0;
}

void bio_set_colnm()
{
    int i;
    if (bio_fmt == BIO_NULL) {
        return;
    } else if (bio_fmt == BIO_HDR) {
        char *p, *q, c;
        for (p = record; *p && isspace(*p); ++p); /* skip leading spaces */
        for (i = 1, q = p; *q; ++q) {
            if (!isspace(*q)) continue;
            c = *q; /* backup the space */
            *q = 0; /* terminate the field */
            set_colnm_aux(p, i);
            *q = c; /* change back */
            ++i;
            for (p = q + 1; *p && isspace(*p); ++p); /* skip contiguous spaces */
            q = p;
        }
        set_colnm_aux(p, i); /* the last column */
    } else {
        for (i = 0; col_defs[bio_fmt][i] != NULL; ++i)
            set_colnm_aux(col_defs[bio_fmt][i], i);
        if (tab_delim[bio_fmt] == 'y') *FS = *OFS = "\t";
    }
}

/**********************
 * Built-in functions *
 **********************/

static char comp_tab[] = {
      0,   1,   2,   3,   4,   5,   6,   7,   8,   9,  10,  11,  12,  13,  14,  15,
     16,  17,  18,  19,  20,  21,  22,  23,  24,  25,  26,  27,  28,  29,  30,  31,
     32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  43,  44,  45,  46,  47,
     48,  49,  50,  51,  52,  53,  54,  55,  56,  57,  58,  59,  60,  61,  62,  63,
     64, 'T', 'V', 'G', 'H', 'E', 'F', 'C', 'D', 'I', 'J', 'M', 'L', 'K', 'N', 'O',
    'P', 'Q', 'Y', 'S', 'A', 'A', 'B', 'W', 'X', 'R', 'Z',  91,  92,  93,  94,  95,
     64, 't', 'v', 'g', 'h', 'e', 'f', 'c', 'd', 'i', 'j', 'm', 'l', 'k', 'n', 'o',
    'p', 'q', 'y', 's', 'a', 'a', 'b', 'w', 'x', 'r', 'z', 123, 124, 125, 126, 127
};

/* The master codon/protein table.
 *
 * http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
 * Tables 7,8,17,18,19,20 are deprecated and do not exist
 * on the NCBI website
 *
 */
const char codon_table[64][25] =
{
/*        1    2    3    4    5    6    7     8     9    10   11   12   13   14   15   16   17    18    19    20   21   22   23   24   25*/
 /*ttt*/{'F', 'F', 'F', 'F', 'F', 'F', '\0', '\0', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', '\0', '\0', '\0', '\0', 'F', 'F', 'F', 'F', 'F'},
 /*ttc*/{'F', 'F', 'F', 'F', 'F', 'F', '\0', '\0', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', '\0', '\0', '\0', '\0', 'F', 'F', 'F', 'F', 'F'},
 /*tta*/{'L', 'L', 'L', 'L', 'L', 'L', '\0', '\0', 'L', 'L', 'L', 'L', 'L', 'L', 'L', 'L', '\0', '\0', '\0', '\0', 'L', 'L', '*', 'L', 'L'},
 /*ttg*/{'L', 'L', 'L', 'L', 'L', 'L', '\0', '\0', 'L', 'L', 'L', 'L', 'L', 'L', 'L', 'L', '\0', '\0', '\0', '\0', 'L', 'L', 'L', 'L', 'L'},
 /*tct*/{'S', 'S', 'S', 'S', 'S', 'S', '\0', '\0', 'S', 'S', 'S', 'S', 'S', 'S', 'S', 'S', '\0', '\0', '\0', '\0', 'S', 'S', 'S', 'S', 'S'},
 /*tcc*/{'S', 'S', 'S', 'S', 'S', 'S', '\0', '\0', 'S', 'S', 'S', 'S', 'S', 'S', 'S', 'S', '\0', '\0', '\0', '\0', 'S', 'S', 'S', 'S', 'S'},
 /*tca*/{'S', 'S', 'S', 'S', 'S', 'S', '\0', '\0', 'S', 'S', 'S', 'S', 'S', 'S', 'S', 'S', '\0', '\0', '\0', '\0', 'S', '*', 'S', 'S', 'S'},
 /*tcg*/{'S', 'S', 'S', 'S', 'S', 'S', '\0', '\0', 'S', 'S', 'S', 'S', 'S', 'S', 'S', 'S', '\0', '\0', '\0', '\0', 'S', 'S', 'S', 'S', 'S'},
 /*tat*/{'Y', 'Y', 'Y', 'Y', 'Y', 'Y', '\0', '\0', 'Y', 'Y', 'Y', 'Y', 'Y', 'Y', 'Y', 'Y', '\0', '\0', '\0', '\0', 'Y', 'Y', 'Y', 'Y', 'Y'},
 /*tac*/{'Y', 'Y', 'Y', 'Y', 'Y', 'Y', '\0', '\0', 'Y', 'Y', 'Y', 'Y', 'Y', 'Y', 'Y', 'Y', '\0', '\0', '\0', '\0', 'Y', 'Y', 'Y', 'Y', 'Y'},
 /*taa*/{'*', '*', '*', '*', '*', 'Q', '\0', '\0', '*', '*', '*', '*', '*', 'Y', '*', '*', '\0', '\0', '\0', '\0', '*', '*', '*', '*', '*'},
 /*tag*/{'*', '*', '*', '*', '*', 'Q', '\0', '\0', '*', '*', '*', '*', '*', '*', 'Q', 'L', '\0', '\0', '\0', '\0', '*', 'L', '*', '*', '*'},
 /*tgt*/{'C', 'C', 'C', 'C', 'C', 'C', '\0', '\0', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', '\0', '\0', '\0', '\0', 'C', 'C', 'C', 'C', 'C'},
 /*tgc*/{'C', 'C', 'C', 'C', 'C', 'C', '\0', '\0', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', '\0', '\0', '\0', '\0', 'C', 'C', 'C', 'C', 'C'},
 /*tga*/{'*', 'W', 'W', 'W', 'W', '*', '\0', '\0', 'W', 'C', '*', '*', 'W', 'W', '*', '*', '\0', '\0', '\0', '\0', 'W', '*', '*', 'W', 'G'},
 /*tgg*/{'W', 'W', 'W', 'W', 'W', 'W', '\0', '\0', 'W', 'W', 'W', 'W', 'W', 'W', 'W', 'W', '\0', '\0', '\0', '\0', 'W', 'W', 'W', 'W', 'W'},
 /*ctt*/{'L', 'L', 'T', 'L', 'L', 'L', '\0', '\0', 'L', 'L', 'L', 'L', 'L', 'L', 'L', 'L', '\0', '\0', '\0', '\0', 'L', 'L', 'L', 'L', 'L'},
 /*ctc*/{'L', 'L', 'T', 'L', 'L', 'L', '\0', '\0', 'L', 'L', 'L', 'L', 'L', 'L', 'L', 'L', '\0', '\0', '\0', '\0', 'L', 'L', 'L', 'L', 'L'},
 /*cta*/{'L', 'L', 'T', 'L', 'L', 'L', '\0', '\0', 'L', 'L', 'L', 'L', 'L', 'L', 'L', 'L', '\0', '\0', '\0', '\0', 'L', 'L', 'L', 'L', 'L'},
 /*ctg*/{'L', 'L', 'T', 'L', 'L', 'L', '\0', '\0', 'L', 'L', 'L', 'S', 'L', 'L', 'L', 'L', '\0', '\0', '\0', '\0', 'L', 'L', 'L', 'L', 'L'},
 /*cct*/{'P', 'P', 'P', 'P', 'P', 'P', '\0', '\0', 'P', 'P', 'P', 'P', 'P', 'P', 'P', 'P', '\0', '\0', '\0', '\0', 'P', 'P', 'P', 'P', 'P'},
 /*ccc*/{'P', 'P', 'P', 'P', 'P', 'P', '\0', '\0', 'P', 'P', 'P', 'P', 'P', 'P', 'P', 'P', '\0', '\0', '\0', '\0', 'P', 'P', 'P', 'P', 'P'},
 /*cca*/{'P', 'P', 'P', 'P', 'P', 'P', '\0', '\0', 'P', 'P', 'P', 'P', 'P', 'P', 'P', 'P', '\0', '\0', '\0', '\0', 'P', 'P', 'P', 'P', 'P'},
 /*ccg*/{'P', 'P', 'P', 'P', 'P', 'P', '\0', '\0', 'P', 'P', 'P', 'P', 'P', 'P', 'P', 'P', '\0', '\0', '\0', '\0', 'P', 'P', 'P', 'P', 'P'},
 /*cat*/{'H', 'H', 'H', 'H', 'H', 'H', '\0', '\0', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', '\0', '\0', '\0', '\0', 'H', 'H', 'H', 'H', 'H'},
 /*cac*/{'H', 'H', 'H', 'H', 'H', 'H', '\0', '\0', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', '\0', '\0', '\0', '\0', 'H', 'H', 'H', 'H', 'H'},
 /*caa*/{'Q', 'Q', 'Q', 'Q', 'Q', 'Q', '\0', '\0', 'Q', 'Q', 'Q', 'Q', 'Q', 'Q', 'Q', 'Q', '\0', '\0', '\0', '\0', 'Q', 'Q', 'Q', 'Q', 'Q'},
 /*cag*/{'Q', 'Q', 'Q', 'Q', 'Q', 'Q', '\0', '\0', 'Q', 'Q', 'Q', 'Q', 'Q', 'Q', 'Q', 'Q', '\0', '\0', '\0', '\0', 'Q', 'Q', 'Q', 'Q', 'Q'},
 /*cgt*/{'R', 'R', 'R', 'R', 'R', 'R', '\0', '\0', 'R', 'R', 'R', 'R', 'R', 'R', 'R', 'R', '\0', '\0', '\0', '\0', 'R', 'R', 'R', 'R', 'R'},
 /*cgc*/{'R', 'R', 'R', 'R', 'R', 'R', '\0', '\0', 'R', 'R', 'R', 'R', 'R', 'R', 'R', 'R', '\0', '\0', '\0', '\0', 'R', 'R', 'R', 'R', 'R'},
 /*cga*/{'R', 'R', 'R', 'R', 'R', 'R', '\0', '\0', 'R', 'R', 'R', 'R', 'R', 'R', 'R', 'R', '\0', '\0', '\0', '\0', 'R', 'R', 'R', 'R', 'R'},
 /*cgg*/{'R', 'R', 'R', 'R', 'R', 'R', '\0', '\0', 'R', 'R', 'R', 'R', 'R', 'R', 'R', 'R', '\0', '\0', '\0', '\0', 'R', 'R', 'R', 'R', 'R'},
 /*att*/{'I', 'I', 'I', 'I', 'I', 'I', '\0', '\0', 'I', 'I', 'I', 'I', 'I', 'I', 'I', 'I', '\0', '\0', '\0', '\0', 'I', 'I', 'I', 'I', 'I'},
 /*atc*/{'I', 'I', 'I', 'I', 'I', 'I', '\0', '\0', 'I', 'I', 'I', 'I', 'I', 'I', 'I', 'I', '\0', '\0', '\0', '\0', 'I', 'I', 'I', 'I', 'I'},
 /*ata*/{'I', 'M', 'M', 'I', 'M', 'I', '\0', '\0', 'I', 'I', 'I', 'I', 'M', 'I', 'I', 'I', '\0', '\0', '\0', '\0', 'M', 'I', 'I', 'I', 'I'},
 /*atg*/{'M', 'M', 'M', 'M', 'M', 'M', '\0', '\0', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', '\0', '\0', '\0', '\0', 'M', 'M', 'M', 'M', 'M'},
 /*act*/{'T', 'T', 'T', 'T', 'T', 'T', '\0', '\0', 'T', 'T', 'T', 'T', 'T', 'T', 'T', 'T', '\0', '\0', '\0', '\0', 'T', 'T', 'T', 'T', 'T'},
 /*acc*/{'T', 'T', 'T', 'T', 'T', 'T', '\0', '\0', 'T', 'T', 'T', 'T', 'T', 'T', 'T', 'T', '\0', '\0', '\0', '\0', 'T', 'T', 'T', 'T', 'T'},
 /*aca*/{'T', 'T', 'T', 'T', 'T', 'T', '\0', '\0', 'T', 'T', 'T', 'T', 'T', 'T', 'T', 'T', '\0', '\0', '\0', '\0', 'T', 'T', 'T', 'T', 'T'},
 /*acg*/{'T', 'T', 'T', 'T', 'T', 'T', '\0', '\0', 'T', 'T', 'T', 'T', 'T', 'T', 'T', 'T', '\0', '\0', '\0', '\0', 'T', 'T', 'T', 'T', 'T'},
 /*aat*/{'N', 'N', 'N', 'N', 'N', 'N', '\0', '\0', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', '\0', '\0', '\0', '\0', 'N', 'N', 'N', 'N', 'N'},
 /*aac*/{'N', 'N', 'N', 'N', 'N', 'N', '\0', '\0', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', '\0', '\0', '\0', '\0', 'N', 'N', 'N', 'N', 'N'},
 /*aaa*/{'K', 'K', 'K', 'K', 'K', 'K', '\0', '\0', 'N', 'K', 'K', 'K', 'K', 'N', 'K', 'K', '\0', '\0', '\0', '\0', 'N', 'K', 'K', 'K', 'K'},
 /*aag*/{'K', 'K', 'K', 'K', 'K', 'K', '\0', '\0', 'K', 'K', 'K', 'K', 'K', 'K', 'K', 'K', '\0', '\0', '\0', '\0', 'K', 'K', 'K', 'K', 'K'},
 /*agt*/{'S', 'S', 'S', 'S', 'S', 'S', '\0', '\0', 'S', 'S', 'S', 'S', 'S', 'S', 'S', 'S', '\0', '\0', '\0', '\0', 'S', 'S', 'S', 'S', 'S'},
 /*agc*/{'S', 'S', 'S', 'S', 'S', 'S', '\0', '\0', 'S', 'S', 'S', 'S', 'S', 'S', 'S', 'S', '\0', '\0', '\0', '\0', 'S', 'S', 'S', 'S', 'S'},
 /*aga*/{'R', '*', 'R', 'R', 'S', 'R', '\0', '\0', 'S', 'R', 'R', 'R', 'G', 'S', 'R', 'R', '\0', '\0', '\0', '\0', 'S', 'R', 'R', 'S', 'R'},
 /*agg*/{'R', '*', 'R', 'R', 'S', 'R', '\0', '\0', 'S', 'R', 'R', 'R', 'G', 'S', 'R', 'R', '\0', '\0', '\0', '\0', 'S', 'R', 'R', 'K', 'R'},
 /*gtt*/{'V', 'V', 'V', 'V', 'V', 'V', '\0', '\0', 'V', 'V', 'V', 'V', 'V', 'V', 'V', 'V', '\0', '\0', '\0', '\0', 'V', 'V', 'V', 'V', 'V'},
 /*gtc*/{'V', 'V', 'V', 'V', 'V', 'V', '\0', '\0', 'V', 'V', 'V', 'V', 'V', 'V', 'V', 'V', '\0', '\0', '\0', '\0', 'V', 'V', 'V', 'V', 'V'},
 /*gta*/{'V', 'V', 'V', 'V', 'V', 'V', '\0', '\0', 'V', 'V', 'V', 'V', 'V', 'V', 'V', 'V', '\0', '\0', '\0', '\0', 'V', 'V', 'V', 'V', 'V'},
 /*gtg*/{'V', 'V', 'V', 'V', 'V', 'V', '\0', '\0', 'V', 'V', 'V', 'V', 'V', 'V', 'V', 'V', '\0', '\0', '\0', '\0', 'V', 'V', 'V', 'V', 'V'},
 /*gct*/{'A', 'A', 'A', 'A', 'A', 'A', '\0', '\0', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', '\0', '\0', '\0', '\0', 'A', 'A', 'A', 'A', 'A'},
 /*gcc*/{'A', 'A', 'A', 'A', 'A', 'A', '\0', '\0', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', '\0', '\0', '\0', '\0', 'A', 'A', 'A', 'A', 'A'},
 /*gca*/{'A', 'A', 'A', 'A', 'A', 'A', '\0', '\0', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', '\0', '\0', '\0', '\0', 'A', 'A', 'A', 'A', 'A'},
 /*gcg*/{'A', 'A', 'A', 'A', 'A', 'A', '\0', '\0', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', '\0', '\0', '\0', '\0', 'A', 'A', 'A', 'A', 'A'},
 /*gat*/{'D', 'D', 'D', 'D', 'D', 'D', '\0', '\0', 'D', 'D', 'D', 'D', 'D', 'D', 'D', 'D', '\0', '\0', '\0', '\0', 'D', 'D', 'D', 'D', 'D'},
 /*gac*/{'D', 'D', 'D', 'D', 'D', 'D', '\0', '\0', 'D', 'D', 'D', 'D', 'D', 'D', 'D', 'D', '\0', '\0', '\0', '\0', 'D', 'D', 'D', 'D', 'D'},
 /*gaa*/{'E', 'E', 'E', 'E', 'E', 'E', '\0', '\0', 'E', 'E', 'E', 'E', 'E', 'E', 'E', 'E', '\0', '\0', '\0', '\0', 'E', 'E', 'E', 'E', 'E'},
 /*gag*/{'E', 'E', 'E', 'E', 'E', 'E', '\0', '\0', 'E', 'E', 'E', 'E', 'E', 'E', 'E', 'E', '\0', '\0', '\0', '\0', 'E', 'E', 'E', 'E', 'E'},
 /*ggt*/{'G', 'G', 'G', 'G', 'G', 'G', '\0', '\0', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', '\0', '\0', '\0', '\0', 'G', 'G', 'G', 'G', 'G'},
 /*ggc*/{'G', 'G', 'G', 'G', 'G', 'G', '\0', '\0', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', '\0', '\0', '\0', '\0', 'G', 'G', 'G', 'G', 'G'},
 /*gga*/{'G', 'G', 'G', 'G', 'G', 'G', '\0', '\0', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', '\0', '\0', '\0', '\0', 'G', 'G', 'G', 'G', 'G'},
 /*ggg*/{'G', 'G', 'G', 'G', 'G', 'G', '\0', '\0', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', '\0', '\0', '\0', '\0', 'G', 'G', 'G', 'G', 'G'},
};

const unsigned char ntval4[256] = {
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 2, 4, 1, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 0, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 2, 4, 1, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 0, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4
    };

char bio_lookup_codon(char *dna, int table)
{
    int ix;
    int i;

    ix = 0;
    for (i=0; i<3; ++i)
    {
        int bv = ntval4[(int)dna[i]];
        if (bv > 3)
            return 'X';
        ix = (ix<<2) + bv;
    }
    return codon_table[ix][table];
}

void bio_translate(char *dna, char *out, int table)
{
    switch(table) {
        case 6:
        case 7:
        case 16:
        case 17:
        case 18:
        case 19:
            WARNING("Translation table %d does not exist. Setting to 1\n", table);
            table = 0;
            break;
        default:
            if(table > 24) {
                WARNING("Translation table %d does not exist. Setting to 1\n", table);
                table = 0;
            }
            break;
    }

    int i;
    int dnaSize;
    int protSize = 0;

    dnaSize = strlen(dna);
    for (i=0; i<dnaSize-2; i+=3)
    {
        out[protSize++] = bio_lookup_codon(dna+i, table);
    }
    for(i=protSize; i < dnaSize; ++i)
        out[i] = '\0';
}

void bio_attribute(Cell * x, Cell * ap, Cell * posp, Cell * y, char kw_delimiter, int del_val_quotes) {
    /* attribute(x, array [, pos_array])
            x, which is a string with the tags or attributes field
              of either gff, or gtf format.
            array, which is the array created by parsing the string.
            optional pos_array, key is field number and value is the field's key in array.

            will return the number of keys in the array in y.

            kw_delimiter '=' for gff ' ' for gtf, del_val_quotes true for gtf
    */

    char *origS, *s, sep, sep2, *sep2_loc, *key, *value, temp;
    int inquote;
    const char QUOTE = '"';
    origS = s = strdup(getsval(x));
    sep = ';'; sep2 = kw_delimiter;
    int n = 0;  // number of fields held in n

    if (*s == '.' && (*(s+1)=='\0' || *(s+1)==' ')) // empty field can be represented by a dot
        *s = '\0';  // drop through to report 0 fields

    while (*s != '\0') {  // make sure not to process empty string and ignore semi-colon at end of string
        if (*s == sep || *s == sep2 || *s == ' ' || *s == '\n') {  // handle doubled semi-colons, prefix spaces and partial handling for other malformed items
            s++;
            continue;
        }
        key = s;

        sep2_loc = NULL; inquote = 0; // set value delimiter loc and handle sep char in quotes
        while ( (*s != sep || inquote) && *s != '\n' && *s != '\0' ) {
            if (*s == sep2 && sep2_loc==NULL)  // first equal sign encountered is value separator
                sep2_loc = s;
            else if (*s == QUOTE)
                inquote = !inquote;
            s++;
        }

        if (sep2_loc == NULL) // if no equal sign, then invalid format
            continue;

        n++;  // increment count of fields

        /*we save the character which should be sep or '\0', then change it
        to the null character to terminate the string. After we set the
        string in the dictionary we set the character back*/
        temp = *s;
        *s = '\0';

        /* ditto with value separator, so key ends with null char */
        *sep2_loc = '\0';

        // replace any spaces before sep2_loc with nulls, so key has ending spaces trimmed
        for (char* pc=(sep2_loc-1); pc > key && *pc==' '; pc--)
            *pc = '\0';

        value = sep2_loc; value++;
        while (*value == ' ') value++;  // skip spaces after equal sign, so value has beginning spaces trimmed

        // replace any spaces before separator with nulls, so value has ending spaces trimmed
        char *valend;
        for (valend=(s-1); valend > value && *valend==' '; valend--)
            *valend = '\0';

        // for gtf files we remove the quotes around the value
        if (del_val_quotes && *value == QUOTE) {
            value++;
            if (valend >= value && *valend == QUOTE)
                *valend = '\0';
        }

        // bio_attribute_inner(t, &key, &value, &temp2);
        if (is_number(value))
            setsymtab(key, value, atof(value), STR|NUM, (Array *) ap->sval);
        else
            setsymtab(key, value, 0.0, STR, (Array *) ap->sval);

        if (posp) {
            char numstr[50];
            snprintf(numstr, sizeof(numstr), "%d", n);
            setsymtab(numstr, key, 0.0, STR, (Array *) posp->sval);
        }

        *sep2_loc = sep2;
        *s = temp;

        if (*s++ == '\0')
            break;
    }
    free(origS);
    if (y != NULL) {  // y is the variable used to return the result, in this case the number of attribute fields
        y->tval = NUM;
        y->fval = n;
        setfval(y, (Awkfloat) n);
    }
}

static float q_int2real[128];

// for treating N and/or YR as wildcards in edit_dist() call
EdlibEqualityPair N_maps[4]   = {{'N', 'A'}, {'N', 'C'}, {'N', 'G'}, {'N', 'T'}};
EdlibEqualityPair YR_maps[4]  = {{'Y', 'C'}, {'Y', 'T'}, {'R', 'A'}, {'R', 'G'}};
EdlibEqualityPair NYR_maps[8] = {{'N', 'A'}, {'N', 'C'}, {'N', 'G'}, {'N', 'T'}, {'Y', 'C'}, {'Y', 'T'}, {'R', 'A'}, {'R', 'G'}};

#define tempfree(x)	  if (istemp(x)) tfree(x); else

Cell *bio_func(int f, Cell *x, Node **a)
{
    Cell *y, *z;
    y = gettemp();
    if (f == BIO_FAND) {
        if (a[1]->nnext == 0) {
            WARNING("and requires two arguments; returning 0.0");
            setfval(y, 0.0);
        } else {
            z = execute(a[1]->nnext);
            setfval(y, (Awkfloat)((long)getfval(x) & (long)getfval(z))); /* FIXME: does (long) always work??? */
            tempfree(z);
        }
    } else if (f == BIO_FOR) {
        if (a[1]->nnext == 0) {
            WARNING("or requires two arguments; returning 0.0");
            setfval(y, 0.0);
        } else {
            z = execute(a[1]->nnext);
            setfval(y, (Awkfloat)((long)getfval(x) | (long)getfval(z)));
            tempfree(z);
        }
    } else if (f == BIO_FXOR) {
        if (a[1]->nnext == 0) {
            WARNING("xor requires two arguments; returning 0.0");
            setfval(y, 0.0);
        } else {
            z = execute(a[1]->nnext);
            setfval(y, (Awkfloat)((long)getfval(x) ^ (long)getfval(z)));
            tempfree(z);
        }
    } else if (f == BIO_FREVERSE) {
        char *buf = getsval(x);
        int i, l, tmp;
        l = strlen(buf);
        for (i = 0; i < l>>1; ++i)
            tmp = buf[i], buf[i] = buf[l-1-i], buf[l-1-i] = tmp;
        setsval(y, buf);
    } else if (f == BIO_FREVCOMP) {
        char *buf;
        int i, l, tmp;
        buf = getsval(x);
        l = strlen(buf);
        for (i = 0; i < l>>1; ++i)
            tmp = comp_tab[(int)buf[i]], buf[i] = comp_tab[(int)buf[l-1-i]], buf[l-1-i] = tmp;
        if (l&1) buf[l>>1] = comp_tab[(int)buf[l>>1]];
        setsval(y, buf);
    } else if (f == BIO_FGC) {
        char *buf;
        int i, l, gc = 0;
        buf = getsval(x);
        l = strlen(buf);
        if (l) { /* don't try for empty strings */
            for (i = 0; i < l; ++i)
                if (buf[i] == 'g' || buf[i] == 'c' ||
                    buf[i] == 'G' || buf[i] == 'C')
                    gc++;
            setfval(y, (Awkfloat)gc / l);
        }
    } else if (f == BIO_FMEANQUAL) {
        char *buf;
        int i, l, total_qual = 0;
        buf = getsval(x);
        l = strlen(buf);
        if (l) { /* don't try for empty strings */
            for (i = 0; i < l; ++i)
                total_qual += buf[i] - 33;
            setfval(y, (Awkfloat)total_qual / l);
        }
    } else if (f == BIO_FTRIMQ) {
        char *buf;
        double thres = 0.05, s, max;
        int i, l, tmp, beg, end;
        Cell *u = 0, *v = 0;
        if (a[1]->nnext) {
            u = execute(a[1]->nnext); /* begin */
            if (a[1]->nnext->nnext) {
                v = execute(a[1]->nnext->nnext); /* end */
                if (a[1]->nnext->nnext->nnext) {
                    z = execute(a[1]->nnext->nnext->nnext);
                    thres = getfval(z); /* user defined threshold */
                    tempfree(z);
                }
            }
        }
        buf = getsval(x);
        l = strlen(buf);
        if (q_int2real[0] == 0.) /* to initialize */
            for (i = 0; i < 128; ++i)
                q_int2real[i] = pow(10., -(i - 33) / 10.);
        for (i = 0, beg = tmp = 0, end = l, s = max = 0.; i < l; ++i) {
            int q = buf[i];
            if (q < 36) q = 36;
            if (q > 127) q = 127;
            s += thres - q_int2real[q];
            if (s > max) max = s, beg = tmp, end = i + 1;
            if (s < 0) s = 0, tmp = i + 1;
        }
        if (u) { setfval(u, beg); tempfree(u); } /* 1-based position; as substr() is 1-based. */
        if (v) { setfval(v, end); tempfree(v); }
    } else if (f == BIO_FQUALCOUNT) {
        if (a[1]->nnext == 0) {
            WARNING("qualcount requires two arguments; returning 0.0");
            setfval(y, 0.0);
        } else {
            char *buf;
            int i, l, thres, cnt = 0;
            buf = getsval(x);
            l = strlen(buf);
            z = execute(a[1]->nnext); /* threshold */
            thres = (int)(getfval(z) + .499);
            for (i = 0; i < l; ++i)
                if (buf[i] - 33 >= thres) ++cnt;
            setfval(y, (Awkfloat)cnt);
        }
    } else if (f == BIO_TRANSLATE) {
        int transtable = 0;
        char *buf;
        char *out;
        if (a[1]->nnext != 0) {
            z = execute(a[1]->nnext);
            transtable = (int) getfval(z);
            /*convert to 0-based indexing*/
            --transtable;
            tempfree(z);
        }
        buf = getsval(x);
        out = calloc(strlen(buf), sizeof(char));
        bio_translate(buf, out, transtable);
        setsval(y, out);
        free(out);
    } else if (f == BIO_GFFATTR || f == BIO_GTFATTR) {
        Cell * ap = NULL; Cell * posp = NULL;
        if (a[1]->nnext != 0) {
            ap = execute(a[1]->nnext);
            if (a[1]->nnext->nnext != 0) {  // 3rd arg names array for holding field positions pos[1], pos[2] etc.
                posp = execute(a[1]->nnext->nnext);
                freesymtab(posp);
                posp->tval &= ~STR; posp->tval |= ARR;
                posp->sval = (char *) makesymtab(NSYMTAB);
            }
        } else {
            char usage[150];
            strcpy(usage, "gffattr(attr_str, array[, pos_array]) requires at least two arguments and allows an optional third");
            if (f == BIO_GTFATTR)
                usage[1] = 't';
            FATAL(usage);
        }
        freesymtab(ap);
        ap->tval &= ~STR;
        ap->tval |= ARR;
        ap->sval = (char *) makesymtab(NSYMTAB);
        char kw_delimiter = (f == BIO_GTFATTR) ? ' ' : '=';
        int del_val_quotes = (f == BIO_GTFATTR);
        bio_attribute(x, ap, posp, y, kw_delimiter, del_val_quotes);
    } else if (f == BIO_FSYSTIME) { /* 12Aug2019 JBH_CAS add systime() that gawk has had for awhile */
        time_t lclock;
        (void) time(& lclock);
        setfval(y, (Awkfloat)lclock);
    } else if (f == BIO_FSETAT) {  /* 06Mar2019 JBH_CAS set characters by position of string (eg seq) */ /* setat(str,pos,vals[,optional_repeat_count])  */
        Cell *u = 0, *v = 0;
        char *buf = getsval(x);
        int l = strlen(buf); int WARN = (l==0);
        if (l) { /* don't try for empty strings */
            if (a[1]->nnext == 0 || a[1]->nnext->nnext == 0) {
                WARN = 1;
            } else {
                u = execute(a[1]->nnext);  /* 1-indexed position */
                v = execute(a[1]->nnext->nnext);  /* string to overlay at position */
                int ix = -1 + (int)getfval(u);
                char *rep = getsval(v);
                int repeat_count = 1;
                if (a[1]->nnext->nnext->nnext) {
                    z = execute(a[1]->nnext->nnext->nnext);
                    repeat_count = (int)getfval(z);
                    tempfree(z);
                }
                if (rep && *rep && ix >= 0 && ix < l) {
                    int i = ix;
                    while (repeat_count-- > 0 && buf[i]) {
                        int r = 0;  /* restart at replacement begin every repeat loop  */
                        for (; buf[i] && rep[r]; i++, r++) {
                            buf[i] = rep[r];
                        }
                    }
                }
            }
        }
        if (WARN) WARNING("setat requires 3 or 4 arguments: str,pos,replacement[, optional repeat_count]");
        if(u!=0){tempfree(u);u=0;} if(v!=0){tempfree(v);v=0;}
        setsval(y, buf);
    } else if (f == BIO_FHAMMING) {  // 26MAR2020 JBH add hamming(pat, text[, text_pos, case_sensitive, N_wildcard ])
        Cell *u = 0, *v = 0;
        int diff = -1, N_wildcard = 0;
        char *pat = getsval(x);
        int compare_len = strlen(pat);
        if (compare_len < 1 || a[1]->nnext == 0) {
            WARNING("hamming requires 2 to 5 arguments: pattern,text [, text_pos: (1_indexed)default 1 [, case_sensitive: true [, N_wildcard: false] ]]");
        } else {
            int case_sensitive = 1, text_ofs = 0;
            u = execute(a[1]->nnext);  /* text string */
            char *text = getsval(u);
            if (a[1]->nnext->nnext != 0) { /* optional text_pos arg, defaults to 1, ie text_ofs 0 */
                v = execute(a[1]->nnext->nnext);  /* 1-indexed text offset */
                text_ofs = -1 + (int)getfval(v);  /* convert awk 1-indexed pos to C 0-indexed offset */
            }
            if (text_ofs < 0) text_ofs = 0;
            text += text_ofs;  /* note: for performance we are not checking bounds which is dangerous in the wild */

            if (a[1]->nnext->nnext && a[1]->nnext->nnext->nnext) { /* if optional case_sensitive arg is there, use it */
                z = execute(a[1]->nnext->nnext->nnext);
                case_sensitive = (0 != (int)getfval(z));
                tempfree(z);
                if (a[1]->nnext->nnext->nnext->nnext) { /* if optional N_wildcard arg is there, use it */
                    z = execute(a[1]->nnext->nnext->nnext->nnext);
                    N_wildcard = (0 != (int)getfval(z)); // non-zero makes N_wildcard
                    tempfree(z);
                }
            }

            diff = 0;
            int to_go = compare_len+1;
            if (N_wildcard==0) {
                if (case_sensitive) {
                    while (--to_go && *text) {
                        diff += *pat != *text;
                        pat++; text++;
                   }
                }
                else {
                   while (--to_go && *text) {
                        diff += toupper(*pat) != toupper(*text);
                        pat++; text++;
                   }
                }
            } else {
                while (--to_go && *text) {
                    char p = *pat, t = *text;
                    if (case_sensitive) {
                        p = toupper(p); t = toupper(t);
                    }
                    if (p!='N' && t!='N')
                        diff += (p != t);
                    pat++; text++;
                }
            }
            diff += to_go;
        }
        if(u!=0){tempfree(u);u=0;} if(v!=0){tempfree(v);v=0;}
        setfval(y, diff);
    } else if (f == BIO_FEDLIB) {  // 26MAR2020 JBH add edit_dist using edlib
        // edit_dist(max_dist, str1, str1_match_len, str2, [str2_len[, mode=EDLIB_MODE_SHW (ie 1)]]

        // can use 0 or -1 for length args and we will do strlen() here
        // EDLIB_MODE_NW  (0) computes edit distance of both strings in entirety (global)
        // EDLIB_MODE_SHW (1) gives distance of str1 to prefix of str2
        // EDLIB_MODE_HW  (2) finds best matches of str1 in str2 (infix or local)
        // if mode value has 10 added to it, we show standard CIGAR, if 20 or more extended CIGAR is shown
        Cell *u = 0, *v = 0, *w = 0;  /* for other args, min 4 up to 6 */
        #define EDBLEN 100
        char edit_dist_buf[EDBLEN] = "-1";
        int mode = EDLIB_MODE_SHW; /* default to prefix mode. works well in concert with a str1_match_len shorter than str1 */
        int task = EDLIB_TASK_LOC; /* gets aligment for cigar when mode arg > 9 EDLIB_TASK_PATH */
        int cigar_type = EDLIB_CIGAR_STANDARD; /* EDLIB_CIGAR_STANDARD if 10<=mode<=19 or EDLIB_CIGAR_EXTENDED if mode >= 20 */
        int N_wildcard = 0, YR_wildcard = 0;

        int max_editdist = (int)getfval(x);  /* -1 means no max set, can be very expensive on long str2 and infix mode */
        int WARN = !(a[1]->nnext && a[1]->nnext->nnext && a[1]->nnext->nnext->nnext); /* args: str1, str1_match_len, str2 */
        if (! WARN) {
            u = execute(a[1]->nnext);  /* str1 */
            char *str1 = getsval(u);
            w = execute(a[1]->nnext->nnext); /* str1_match_len */
            int slen1 = (int)getfval(w); tempfree(w); w=0;
            if(slen1 < 1) /* if 0 or -1 passed we'll figure it out here */
               slen1 = strlen(str1);

            v = execute(a[1]->nnext->nnext->nnext); /* str2 */
            char* str2 = getsval(v);
            int slen2 = 0;

            if (a[1]->nnext->nnext->nnext->nnext) { /* optional 5th arg: slen2 */
                w = execute(a[1]->nnext->nnext->nnext->nnext);
                slen2 = (int)getfval(w); tempfree(w); w=0;

                if (a[1]->nnext->nnext->nnext->nnext->nnext) { /* optional 6th arg: mode */
                    w = execute(a[1]->nnext->nnext->nnext->nnext->nnext);
                    int mval = (int)getfval(w); tempfree(w); w=0;
                    if (mval >= 10) { /* mode + 10 means show CIGAR for alignment */
                        if (mval >= 20) cigar_type = EDLIB_CIGAR_EXTENDED;
                        mval = mval % 10;
                        task = EDLIB_TASK_PATH;
                    }
                    if (mval != EDLIB_MODE_SHW) { /* 0 is complete match NW, 1 is prefix match SHW, any other is infix match HW */
                        mode = (mval==EDLIB_MODE_NW) ? EDLIB_MODE_NW : EDLIB_MODE_HW;
                    }
                    if (a[1]->nnext->nnext->nnext->nnext->nnext->nnext) { /* optional 7th arg: equality pair types */
                        w = execute(a[1]->nnext->nnext->nnext->nnext->nnext->nnext);
                        int equalities_flag = (int)getfval(w); tempfree(w); w=0;
                        N_wildcard = equalities_flag & 1; // true if 1 bit set
                        YR_wildcard = equalities_flag & 2; // true if 2 bit set
                    }
                }
            }
            if (slen2 < 1) { /* get slen2 based on actual length; no arg 5 or trigger this by passing in 0 or -1 as str2_len arg */
                slen2 = strlen(str2);
            }
            if (slen1 < 1 || slen2 < 1) {
                char* msg = (slen1 < 1) ? "str1_match_len must be greater than or equal to 1" : "str2 empty";
                WARNING(msg);
            }
            else {
                EdlibEqualityPair *addtlEqualities = NULL; int equalities_len = 0;
                if (N_wildcard || YR_wildcard) {
                    addtlEqualities = (N_wildcard && YR_wildcard) ? NYR_maps : (N_wildcard) ? N_maps : YR_maps;
                    equalities_len  = (N_wildcard && YR_wildcard) ? 8 : 4;
                }

                EdlibAlignConfig edlibConfig = edlibNewAlignConfig(max_editdist, mode, task, addtlEqualities, equalities_len);

                EdlibAlignResult result = edlibAlign(str1, slen1, str2, slen2, edlibConfig);
                if (result.status == EDLIB_STATUS_OK) {
                    char temp[50];
                    int start = -1, end = -1;
                    sprintf(edit_dist_buf, "%d", result.editDistance);
                    if (result.alignment) {
                        char* cigar = edlibAlignmentToCigar(result.alignment, result.alignmentLength, cigar_type);
                        sprintf(temp, " %s", cigar);
                        strcat(edit_dist_buf, temp);
                        free(cigar);
                    }
                    for (int i=0; i < result.numLocations; i++) { // 1 index start and end locations for awk ouput
                        start = result.startLocations[i] + 1;
                        end = result.endLocations[i] + 1;
                        sprintf(temp, " %d %d", start, end);
                        if(strlen(edit_dist_buf)+strlen(temp)+1 < EDBLEN)
                            strcat(edit_dist_buf, temp);
                    }
                }
                edlibFreeAlignResult(result);
            }
        }
        if (WARN) {
            WARNING("edit_dist requires 4 to 7 arguments: max_editdist, str1, str1_match_len, str2[, str2_len [, mode: default 1 [, flags]]]\n"
                    "                  mode: 0 complete match, 1 prefix match, 2 infix match (add 10 or 20 for CIGAR). Can use string len -1 for full length.\n"
                    "                  flags: 1 N matches ACTG, 2 Y matches CT, R matches AG, 3 both.");
        }

        if(u!=0){tempfree(u);u=0;} if(v!=0){tempfree(v);v=0;} if(w!=0){tempfree(w);w=0;}
        setsval(y, edit_dist_buf); /* return string with edit_distance start_loc end_loc */
    }
    else if (f == BIO_ADAPATEND) {
        Cell *u = 0;
        char* seq_to_chk = getsval(x);
        char* adapter = 0;
        char match_info_buf[60] = {'0'};
        int WARN = 0;
        if (seq_to_chk==NULL || *seq_to_chk=='\0') {  // setting adapter when 1st arg is empty, eg end_adapter_pos("", "GATCGGAAGAGCACAC")
            if (a[1]->nnext) { // adapter is in 2nd arg, we only look at it when first arg is empty
                u = execute(a[1]->nnext);
                adapter = getsval(u);
                if (adapter != NULL && *adapter != '\0') {
                    free_g_adap_info();
                    g_adap_info = make_adapter_prefix_encodings(adapter);
                    sprintf(match_info_buf, "%d", g_adap_info.adapter_length);  // length of adapter used for prefix check, max 16
                } else { WARN = 1; }
            } else { WARN = 1; }
        }
        else if (g_adap_info.prefix_set != NULL && g_adap_info.adapter_length >= 4) // calling to check seq, eg end_adapter_pos(seq)
        {
            struct readEndMatch match_inf = check_read_end_for_adapter_prefix(seq_to_chk);
            sprintf(match_info_buf, "%d %d %d", match_inf.match_pos, match_inf.len, match_inf.errs);
        } else { WARN = 1; }

        if (WARN) {
            WARNING("end_adapter_pos(\"\", adapter) to set adapter. end_adapter_pos(seq) to check seq suffix against adapter prefix.\n"
                "                  To set adapter call with empty seq, subsequent calls use seq as only argument.\n"
                "                  Returns string with 3 numbers: position of match, len, mismatches (-1 for none)");
        }

        setsval(y, match_info_buf);

    } else if (f == BIO_FMD5) { /* 26May2020 JBH_CAS add md5() to return md5 string for input parm1 */

        char* seq_to_chk = getsval(x);
        size_t len = strlen(seq_to_chk);

        if (len > 0) {
            char* md5_rslt = md5str((unsigned char*)seq_to_chk, len);
            setsval(y, md5_rslt);
            free(md5_rslt); *md5_rslt = '\0';
        } else {
            WARNING("md5 takes a string argument and returns its md5 value.");
        }

    } else if (f == BIO_CHARCOUNT) { /* charcount(str, ar_chars) -- returns val e.g. ar_chars["A"]=173 JBH 21Jul2020 */

        int char_count = 0, counts[256] = {0};
        char num_str[50], char_as_str[2] = {0,0}; /* char_as_str[0] filled with char */
        Cell *ap = 0;

        int WARN = (a[1]->nnext == 0);
        if ( ! WARN ) {
            ap = execute(a[1]->nnext);     /* array name */
            freesymtab(ap);
            ap->tval &= ~STR;
            ap->tval |= ARR;
            ap->sval = (char *) makesymtab(NSYMTAB);

            char* seq = getsval(x);
            if (seq && *seq) { /* we have chars to look at */

                for (char* pchr = seq; *pchr; pchr++) { /* count chars in string, store in our local array */
                   int ch = (int)*pchr;
                   counts[ch] += 1;
                   char_count += (counts[ch]==1);
                }

                /* for those chars with non-zero count, add them into array (aka symtab) */
                int count;
                for (int c=0; c<256; c++) {
                    if ( (count=counts[c]) ) {
                        char_as_str[0] = (char)c;
                        sprintf(num_str, "%d", count);
                        setsymtab(char_as_str, num_str, count, STR|NUM, (Array *) ap->sval);
                    }
                }
            }
        }

        if (WARN)
            WARNING("charcount(str, arr) fills arr with count of each character in str. returns number of different chars.");

        setfval(y, (Awkfloat)char_count);

    } else if (f == BIO_MODSTR) { /* modstr(str, start, length, mod_type) -- mod_type 0 to lowercase, 1 to uppercase, might add ability to use function for mod_type later 27Jul2020 */
        int WARN = !(a[1]->nnext && a[1]->nnext->nnext); /* 2nd, 3rd args: start, length */

        if (!WARN) {
            char *str = getsval(x);
            int slen = 0; /* optional 5th arg lets us set this ourselves, otherwise use strlen() which can be very slow if we call this thousands of times */

            Cell *u = execute(a[1]->nnext); /* start pos */
            Cell *v = execute(a[1]->nnext->nnext); /* length pos */
            int start = (int) getfval(u);
            tempfree(u);
            u = 0;
            int length = (int) getfval(v);
            tempfree(v);
            v = 0;
            int mod_type = 0; /* 0 means lowercase and that is our default type of modification*/
            if (a[1]->nnext->nnext->nnext) {
                Cell *w = execute(a[1]->nnext->nnext->nnext);
                mod_type = (int) getfval(w);
                tempfree(w);
                w = 0;
                if (a[1]->nnext->nnext->nnext->nnext) { // optional 5th arg to pass in length
                    Cell *w = execute(a[1]->nnext->nnext->nnext->nnext);
                    slen = (int) getfval(w);
                    tempfree(w);
                    w = 0;
                }
            }
            if (slen < 1)
                slen = strlen(str);
            //char msg[100]; sprintf(msg, "%s,%d,%d,%d", str, start, length, mod_type); WARNING(msg);

            int ix = start - 1;
            if (ix > -1 && ix < slen) {
                for (char *p = &str[ix]; *p && length--; p++)
                    *p = (mod_type == 0) ? tolower(*p) : toupper(*p);
            }
            setfval(y, (Awkfloat) 0);
        }

        if (WARN)
            WARNING("modstr(str, start, length, [mod_type, str_length]) mod_type 0 to lowercase, 1 to uppercase (default 0). optional str_length faster for multiple calls.");

    } else if (f == BIO_APPLYCHARS) { /* loop over each char of str arg1 and execute 2nd arg for each setting CHAR and ORD variables */
        Cell *pcellchar = setsymtab("CHAR", "", 0.0, STR, symtab); /* setsymtab either creates and returns ptr to var as a Cell, or returns existing one  */
        Cell *pcellord = setsymtab("ORD", "", 0.0, NUM, symtab); /* ordinal value of CHAR */

        int OK = a[1]->nnext && pcellchar && pcellord;
        if (OK) {
            char *str = getsval(x); /* first arg is the string from which the chars are gotten */
            setsval(pcellchar, "A"); /* allocate memory for char_as_str in Cell */
            for (char *pchr = str; *pchr; pchr++) {
                pcellchar->sval[0] = *pchr; // was setsval(pcellchar, char_as_str); which malloc's memory for each CHAR
                setfval(pcellord, (Awkfloat) *pchr);
                for (Node *eval_stmt = a[1]->nnext; eval_stmt; eval_stmt = eval_stmt->nnext) { /* execute over each char using var CHAR */
                    execute(eval_stmt);
                }
            }
        } else
            WARNING("applytochars(str, stmt_or_func). 2nd arg called for each char in str with CHAR and ORD variables set.");

    } /* else: never happens */
    return y;
}

/************************
 * getrec() replacement *
 ************************/

#include <zlib.h> /* FIXME: it would be better to drop this dependency... */
#include "kseq.h"
KSEQ_INIT2(, gzFile, gzread)

static gzFile g_fp;
static kseq_t *g_kseq;
static int g_firsttime = 1, g_is_stdin = 0;
static kstring_t g_str;

int bio_getrec(char **pbuf, int *psize, int isrecord)
{
    extern Awkfloat *ARGC;
    extern int argno, recsize;
    extern char *file;
    extern Cell **fldtab;

    int i, c, saveb0, dret, bufsize = *psize, savesize = *psize;
    char *p, *buf = *pbuf;
    if (g_firsttime) { /* mimicing initgetrec() in lib.c */
        g_firsttime = 0;
        for (i = 1; i < *ARGC; i++) {
            p = getargv(i); /* find 1st real filename */
            if (p == NULL || *p == '\0') {	/* deleted or zapped */
                argno++;
                continue;
            }
            if (!isclvar(p)) {
                setsval(lookup("FILENAME", symtab), p);
                goto getrec_start;
            }
            setclvar(p);	/* a commandline assignment before filename */
            argno++;
        }
        g_fp = gzdopen(fileno(stdin), "r"); /* no filenames, so use stdin */
        g_kseq = kseq_init(g_fp);
        g_is_stdin = 1;
    }

getrec_start:
    if (isrecord) {
        donefld = 0; /* these are defined in lib.c */
        donerec = 1;
    }
    saveb0 = buf[0];
    buf[0] = 0; /* this is effective at the end of file */
    while (argno < *ARGC || g_is_stdin) {
        if (g_kseq == 0) { /* have to open a new file */
            file = getargv(argno);
            if (file == NULL || *file == '\0') { /* deleted or zapped */
                argno++;
                continue;
            }
            if (isclvar(file)) {	/* a var=value arg */
                setclvar(file);
                argno++;
                continue;
            }
            *FILENAME = file;
            if (*file == '-' && *(file+1) == '\0') {
                g_fp = gzdopen(fileno(stdin), "r");
                g_kseq = kseq_init(g_fp);
                g_is_stdin = 1;
            } else {
                if ((g_fp = gzopen(file, "r")) == NULL)
                    FATAL("can't open file %s", file);
                g_kseq = kseq_init(g_fp);
                g_is_stdin = 0;
            }
            setfval(fnrloc, 0.0);
        }
        if (bio_fmt != BIO_FASTX) {
            c = ks_getuntil(g_kseq->f, **RS, &g_str, &dret);
        } else {
            c = kseq_read(g_kseq);
            if (c >= 0) {
                g_str.l = 0;
                g_str.m = g_kseq->name.l + g_kseq->comment.l + g_kseq->seq.l + g_kseq->qual.l + 4;
                kroundup32(g_str.m);
                g_str.s = (char*)realloc(g_str.s, g_str.m);
                for (i = 0; i < g_kseq->name.l; ++i)
                    g_str.s[g_str.l++] = g_kseq->name.s[i];
                g_str.s[g_str.l++] = '\t';
                for (i = 0; i < g_kseq->seq.l; ++i)
                    g_str.s[g_str.l++] = g_kseq->seq.s[i];
                g_str.s[g_str.l++] = '\t';
                for (i = 0; i < g_kseq->qual.l; ++i)
                    g_str.s[g_str.l++] = g_kseq->qual.s[i];
                g_str.s[g_str.l++] = '\t';
                for (i = 0; i < g_kseq->comment.l; ++i)
                    g_str.s[g_str.l++] = g_kseq->comment.s[i];
                g_str.s[g_str.l++] = '\0';
            } else {
                g_str.l = 0;
                if (g_str.s) g_str.s[0] = '\0';
            }
        }
        adjbuf(&buf, &bufsize, g_str.l + 1, recsize, 0, "bio_getrec");
        memcpy(buf, g_str.s, g_str.l + 1);
        if (c >= 0) {	/* normal record */
            if (isrecord) {
                if (freeable(fldtab[0]))
                    xfree(fldtab[0]->sval);
                fldtab[0]->sval = buf;	/* buf == record */
                fldtab[0]->tval = REC | STR | DONTFREE;
                if (is_number(fldtab[0]->sval)) {
                    fldtab[0]->fval = atof(fldtab[0]->sval);
                    fldtab[0]->tval |= NUM;
                }
            }
            setfval(nrloc, nrloc->fval+1);
            setfval(fnrloc, fnrloc->fval+1);
            *pbuf = buf;
            *psize = bufsize;
            return 1;
        }
        /* EOF arrived on this file; set up next */
        if (!g_is_stdin) {
            kseq_destroy(g_kseq);
            gzclose(g_fp);
        }
        g_fp = 0; g_kseq = 0; g_is_stdin = 0;
        argno++;
    }
    buf[0] = saveb0;
    *pbuf = buf;
    *psize = savesize;
    return 0;	/* true end of file */
}