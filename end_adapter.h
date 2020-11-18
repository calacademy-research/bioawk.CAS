// end_adapter is the code for including this function in bioawk_cas
// gcc suffix_match_snippet.c -o suffix_match_snippet -O2 -mavx2 -pedantic-errors -lz

// -std=gnu17 is the default for C code
// you can use -S to generate the assembly equivalent of the code in a .s file (when not using -o)
// -mavx2 gets popcnt, don't need -mpopcnt, that was for older architectures
// -lz links to the zlib, -lm to the math and -lrt to the realtime extensions library

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>

// prep for 128 bit version that can encode 32 bases
static inline unsigned popcnt_u128 (__uint128_t n)
{
    return __builtin_popcountll((uint64_t) (n >> 64))  // bits set in high-order 64bits of 128bit uint
         + __builtin_popcountll((uint64_t) n);         // bits set in low-order 64bits of 128bit uint
}


// Aa == 1 Cc==2 Gg==4 Tt==8 Nn==15 else 0, so 0001 A, 0010 C, 0100 G, 1000 T, 1111 N
// (N is 15 so it plays part as a wildcard when AND'd with A, C, G, or T)
// if pat to compare with has no Ns and seq to compare may have them, then
// hamming dist is given by compare_len - pop_count( pat_4bit & seq_4bit )

const unsigned char nt4bit[256] = {  // each of the values 0, 1, 2, 8, 15 fits in 4 bits aka nybble
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, //  0 - 15
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 16 - 31
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 32 - 47
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 48 - 63
    0, 1, 0, 2, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 15,0, // 64 - 79
    0, 0, 0, 0, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 80 - 95
    0, 1, 0, 2, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 15,0, // 96 - 111
    0, 0, 0, 0, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 112 - 127
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 128 - 143
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 144 - 159
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 160 - 175
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 176 - 191
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 192 - 207
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 208 - 223
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 224 - 239
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0  // 240 - 255
};

const char from4bit[16] = { 'x', 'A', 'C', 0, 'G', 0, 0, 0, 'T', 0, 0, 0, 0, 0, 0, 'N' };

struct adapterPrefixInfo {
    uint64_t *prefix_set;
    int set_size;
    int adapter_length; // 1 - 16, really 4 - 16
} g_adap_info = {NULL, 0, 0};

void free_g_adap_info() {
    if(g_adap_info.prefix_set){
        free(g_adap_info.prefix_set);
        g_adap_info = (struct adapterPrefixInfo){.prefix_set = NULL, .set_size = 0, .adapter_length = 0};
    }
}

// convert an up to 16 nt sequence into a 64 bit uint with 4bits per ACGTN
uint64_t shortseq_to_4bit(const char *shortseq) {
    uint64_t encoded = 0;
    int mx = 16;
    while (*shortseq && mx--) {
        encoded = (encoded << 4) | nt4bit[ (unsigned) *shortseq ];
        shortseq++;
    }
    while (mx--) // we want nt bits to start at leftmost bits of int
        encoded = (encoded << 4);

    return encoded;
}

char *decode_4bit(uint64_t encoded) {
    char *shortseq = (char*)calloc(17, sizeof(char));
    char *cp = shortseq;
    for (int nyb_end=60; nyb_end >= 0; nyb_end -= 4) {
        unsigned i = (encoded>>nyb_end) & 0xF;
        *cp++ = from4bit[i];
    }
    // so we can leave non-ACTGN chars inside the seq (shown as 'x') above has also filled out with x's to end
    // this will trim those x's at end off
    for(cp=(shortseq + strlen(shortseq)-1); cp >= shortseq && *cp==from4bit[0]; cp--) {
        *cp = '\0';
    }
    return shortseq;
}

// we want to check the end of a read for the beginning of an adapter,
// so we check the last 16 nt of the read against the first 16 nt of the adapter seq,
// then the last 15 nt of the read for the first 15 nt of the adapter
// and so-on until the last 4 nt of the read is checked with first 4 adapter nt.
// we stop when an attempt has an acceptable hamming distance.

// to do this we make an array of upto 13 uint64_t's with the 4bit encodings
// of the adapter 16nt prefix, 15nt prefix, etc. we reuse this set for each read...
struct adapterPrefixInfo make_adapter_prefix_encodings(const char* adapter_prefix) {
    struct adapterPrefixInfo result = {NULL, 0, 0};
    
    const int MAXCOMPARE = 4;
    char adap[17] = {0};
    strncpy(adap, adapter_prefix, 16); // so we have a copy to truncate
    
    int len = strlen(adap);
    int iter = 1 + len - MAXCOMPARE;
    if (iter < 1) return result;
    
    result.adapter_length = len;
    uint64_t *prefix_set = (uint64_t*)calloc(iter+1, sizeof(uint64_t));
    for(int i=0; i<iter; i++) {
        prefix_set[i] = shortseq_to_4bit(adap);
        adap[ --len ] = '\0';
    }
    
    result.prefix_set = prefix_set;
    result.set_size = iter;   
    return result;
}

//int error_threshold [17] = {0, 0, 0, 0, 0, 1, 1, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4};
  int match_threshold [17] = {0, 1, 2, 3, 4, 4, 5, 5, 6, 6, 7, 8, 9, 9,10,11,12};
//                            -  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16

// ...then we take the adapter prefix array and a read, encode the
// last 16nt of the read and compute hamming distance for each
// adapter prefix and read suffix, shortening as we go.
// to get mismatch: AND encoded seqs' bits, get popcount and subtract from seq len
struct readEndMatch { int match_pos; int len; int errs; };
struct readEndMatch check_read_end_for_adapter_prefix(const char* read) {
    static int recs = 0;
    struct readEndMatch match_inf = { -1, -1, -1 };
    if (g_adap_info.prefix_set==NULL || g_adap_info.adapter_length < 4)
        return match_inf; // adapter prefix set not initialized, call make_adapter_prefix_encodings() first
    
    recs++;
    int readofs = strlen(read) - g_adap_info.adapter_length;
    if (readofs < 0) return match_inf;

    uint64_t readtail = shortseq_to_4bit(read+readofs);

    // loop through from longest prefix to shortest and check the number of matches
    int matches = 0;
    int *enough = &match_threshold[16];
    uint64_t *prefixes = g_adap_info.prefix_set;
    int curlen = g_adap_info.adapter_length;
    for(int i=0; i < g_adap_info.set_size; i++, curlen--, enough--) {
        uint64_t prf = prefixes[i];
        uint64_t match_bits = prf & readtail;
        matches = __builtin_popcountll(match_bits);
        
        if (matches >= *enough) {
            match_inf.match_pos = readofs+i;
            match_inf.len = curlen;
            match_inf.errs = curlen-matches;            
            break;
        }
        
        readtail = readtail << 4;   // eg, encoded ATGATC becomes TGATC (0x184182 becomes 0x84182)
    }
    return match_inf;
}
