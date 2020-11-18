#!/bin/bash
./bioawk -c fastx 'BEGIN{end_adapter_pos("", "GATCGGAAGAGCACAC")}
   {rslt=end_adapter_pos($seq)}
   int(rslt)>0{print NR": "rslt; r++}
 ' ~/turtle/Tg_CKDL200149914-1a_H7N5LCCX2_L4_R1.fq.gz
