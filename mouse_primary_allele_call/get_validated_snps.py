"""
Script that analyzes snpstatistics files of different genotypes to identify a list of validated snps
"""

import os, sys

nts = {'A':4, 'C':5, 'G':6, 'T':7, 'N':8 }

fhout = open('validated_cast_c57_snps.txt','w')
fhoutX = open('nonsupported_cast_c57_snps.txt','w')

with open('snpstatistics_cast_c57.txt') as cast_c57, open('snpstatistics_cast.txt') as cast,  open('snpstatistics_c57_cast.txt') as c57_cast, open('snpstatistics_c57.txt') as c57:
    for cb, c, bc, b in zip(cast_c57, cast, c57_cast, c57):
        if cb[0] == '#': continue
        cb_parts = cb.strip().split('\t')
        c_parts  = c.strip().split('\t')
        bc_parts = bc.strip().split('\t')
        b_parts  = b.strip().split('\t')
        
        c57_nt = cb_parts[2]
        cast_nt = cb_parts[3]

        c57_pos = []
        for c57n in c57_nt.split(","):
            c57_pos.append(nts[c57n])
        cast_pos= []
        for castn in cast_nt.split(","):
            cast_pos.append(nts[castn])
        
        # criteria 
        # > 90% correct genotype in pure cells
        # or 
        # at least 10% of each genotype in mixed cells
        
        # evaluate pure c57
        if sum(map(float,b_parts[4:8])) > 0:
            c57_reads_in_c57 = sum( float(b_parts[c57p]) for c57p in c57_pos) / sum(map(float,b_parts[4:8]))
        else:
            c57_reads_in_c57 = -1

        # evaluate pure cast
        if sum(map(float,c_parts[4:8])) > 0:
            cast_reads_in_cast = sum( float(c_parts[castp]) for castp in cast_pos) / sum(map(float,c_parts[4:8]))
        else:
            cast_reads_in_cast = -1

        # evaluate alleles in mixed genotype cells
        if  sum(map(float,cb_parts[4:8])) > 0:
            #print (cast_pos)
            cast_reads_in_mixed_cb = sum( float(cb_parts[castp]) for castp in cast_pos) / sum(map(float,cb_parts[4:8]))
            c57_reads_in_mixed_cb = sum( float(cb_parts[c57p]) for c57p in c57_pos) / sum(map(float,cb_parts[4:8]))
        else:
            cast_reads_in_mixed_cb, c57_reads_in_mixed_cb = -1, -1
        if sum(map(float,bc_parts[4:8])) > 0:
            cast_reads_in_mixed_bc = sum( float(bc_parts[castp]) for castp in cast_pos) / sum(map(float,bc_parts[4:8]))
            c57_reads_in_mixed_bc = sum( float(bc_parts[c57p]) for c57p in c57_pos) / sum(map(float,bc_parts[4:8]))
        else:
            cast_reads_in_mixed_bc, c57_reads_in_mixed_bc = -1, -1


        flagged = 1
        if (c57_reads_in_c57 > 0.9 and cast_reads_in_cast > 0.9)  or \
                min(cast_reads_in_mixed_cb, c57_reads_in_mixed_cb) > 0.1 or \
                min(cast_reads_in_mixed_bc, c57_reads_in_mixed_bc) > 0.1:
            flagged = 0

        if flagged:
            fhoutX.write("%s\t%s\t%s\n" % ("\t".join(c_parts[:4]),
                                          flagged,
                                          "\t".join(map(str, [c57_reads_in_c57,
                                                              cast_reads_in_cast,
                                                              cast_reads_in_mixed_cb,
                                                              c57_reads_in_mixed_cb,
                                                              cast_reads_in_mixed_bc,
                                                              c57_reads_in_mixed_bc]))))


        else:
            fhout.write("%s\t%s\t%s\n" % ("\t".join(c_parts[:4]),
                                          flagged,
                                          "\t".join(["%.2f"%val for val in [c57_reads_in_c57,
                                                                            cast_reads_in_cast,
                                                                            cast_reads_in_mixed_cb,
                                                                            c57_reads_in_mixed_cb,
                                                                            cast_reads_in_mixed_bc,
                                                                            c57_reads_in_mixed_bc]])))
            

                               
