#!/usr/bin/env python
# *- coding: utf-8 -*-
# vim: set expandtab tabstop=4 shiftwidth=4 softtabstop=4 textwidth=79:

"""[License: GNU General Public License v3 (GPLv3)]

    fibronectin-splice-variant-detector: counts Fibronectin (FN1) alt. splicing in BAM files
    Copyright (C) 2024  Youri Hoogstrate, Tobias Weiss and Pim French

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.


    You can contact me via the github repository at the following url:
    <https://github.com/yhoogstrate/fibronectin-splice-var-determiner>

    You can e-mail me via 'y.hoogstrate' at the following webmail domain:
    gmail dot com
"""

import sys
import pysam


# search window 7B89
# chr2:215,390,931-215,395,491

# search window 11Q12
# chr2:215,378,835-215,382,551

# SV 24->26 read: A00379:269:HGTK3DSXY:3:2208:30255:24330
# SV 24->25 read: A00379:263:HGNF2DSXY:4:2229:14552:2300





fn1_exons = { # fibronectin
  # 'hg19': { },
  'hg38': {
    "24": ['chr2', 215394528, 215394719, "-"], # ex-24 
    "25": ['chr2', 215392931, 215393203, "-"], # ex-25: 7B89 & strand = neg
    "26": ['chr2', 215391632, 215391814, "-"], # ex-26
 
    "31": ["chr2", 215382212, 215382325, "-"], # ex-31
    "32": ["chr2", 215380811, 215381080, "-"], # ex-32: 11A12
    "33": ["chr2", 215379130, 215379317, "-"]  # ex-33
    }
}


def check_or_update_chr(chrlist, pysam_handle):
    for _ in chrlist:
        _chr = chrlist[_][0]
        if _chr not in pysam_handle.references:
            chrlist[_][0] = chrlist[_][0].replace('chr','')

    return chrlist


def get_splice_junction_positions(alignedsegment):
    """
    https://sourceforge.net/p/samtools/mailman/message/29373646/

    M    BAM_CMATCH    0
    I    BAM_CINS    1
    D    BAM_CDEL    2
    N    BAM_CREF_SKIP    3
    S    BAM_CSOFT_CLIP    4
    H    BAM_CHARD_CLIP    5
    P    BAM_CPAD    6
    =    BAM_CEQUAL    7
    X    BAM_CDIFF    8
    B    BAM_CBACK    9
    """
    out = []
    #out = [None, None]# start, end
    offset = alignedsegment.reference_start
    for cigar in alignedsegment.cigartuples:
        s = [offset,cigar]
        
        # Insertion adds a gap, does not pad
        # soft and hardclip don't add padding
        #               M  D  =  X
        if cigar[0] in [0, 2, 7, 8]:# Match, Insertion, padding, Read match?, mismatch
            offset += cigar[1]
            s.append(offset)
        elif cigar[0] == 9:
            offset -= cigar[1]
        elif cigar[0] == 3: # splice junction starts
            #out[0] = offset
            #out[1] = offset + cigar[1] + 1

            out.append([offset, offset + cigar[1] + 1])
            
            offset += cigar[1]
    
    return out


def extract_FN1_7B89_sv_reads(bam, exons, include_interchromosomal, include_duplicates):
    #t = {'24': {'A->B->C', 'A->B', 'A->C'}, '25': {'A->B->C', 'A->B', 'B->C'}, '26': {'A->B->C', 'A->C', 'B->C'}}
    readnames = {'24': set([]), '25': set([]), '26': set([])}

    fh = pysam.AlignmentFile(bam, "rb")
    exons = check_or_update_chr(exons, fh)

    for exon in exons:
        for read in fh.fetch(exons[exon][0], exons[exon][1], exons[exon][2]):
            if read.get_overlap(exons[exon][1], exons[exon][2]):
                if include_interchromosomal or (not read.is_paired or (read.is_paired and read.next_reference_name in list(set([_[0] for _ in exons.values()])))):
                    if include_duplicates or (not read.is_duplicate):
                        if exon in readnames:
                            readnames[exon].add(read.query_name)

    wt = readnames['24'].intersection(readnames['26']).difference(readnames['25'])
    splice_right_intron = readnames['24'].intersection(readnames['25']).difference(readnames['26'])
    splice_left_intron = readnames['25'].intersection(readnames['26']).difference(readnames['24'])
    splice_both_introns = readnames['24'].intersection(readnames['25'], readnames['26'])
    

    return {'7B89: wt read [24 -> 26]': wt,
            '7B89: sv read [24 -> 25 -> 26]': splice_both_introns,
            '7B89: sv read [24 -> 25]': splice_right_intron,
            '7B89: sv read [25 -> 26]': splice_left_intron,
            '7B89: discrepant read': set([])} # this does not test for discrepancies


def extract_FN1_11A12_sv_reads(bam, exons, include_interchromosomal, include_duplicates):
    readnames = {'31': set([]), '32': set([]), '33': set([])}

    fh = pysam.AlignmentFile(bam, "rb")
    exons = check_or_update_chr(exons, fh)

    for exon in exons:
        for read in fh.fetch(exons[exon][0], exons[exon][1], exons[exon][2]):
            if read.get_overlap(exons[exon][1], exons[exon][2]):
                if include_interchromosomal or (not read.is_paired or (read.is_paired and read.next_reference_name in list(set([_[0] for _ in exons.values()])))):
                    if include_duplicates or (not read.is_duplicate):
                        if exon in readnames:
                            readnames[exon].add(read.query_name)

    wt = readnames['31'].intersection(readnames['33']).difference(readnames['32'])
    splice_right_intron = readnames['31'].intersection(readnames['32']).difference(readnames['33'])
    splice_left_intron = readnames['32'].intersection(readnames['33']).difference(readnames['31'])
    splice_both_introns = readnames['31'].intersection(readnames['32'], readnames['33'])
    

    return {'11A12: wt read [31 -> 33]': wt,
            '11A12: sv read [31 -> 32 -> 33]': splice_both_introns,
            '11A12: sv read [31 -> 32]': splice_right_intron,
            '11A12: sv read [32 -> 33]': splice_left_intron,
            '11A12: discrepant read': set([])} # this does not test for discrepancies


def extract_FN1_7B89_sv_reads_based_on_sjs(bam, exons, include_interchromosomal, include_duplicates):
    set_wt = set()# readnames of those that splice from exon 24 to 26
    set_sv1 = set()# readnames of those that splice from exon 24 to 25
    set_sv2 = set()# readnames of those that splice from exon 25 to 26
    
    fh = pysam.AlignmentFile(bam, "rb")
    exons = check_or_update_chr(exons, fh)

    exon = "24"
    for read in fh.fetch(exons[exon][0], exons[exon][1], exons[exon][2]):
        if read.get_overlap(exons[exon][1], exons[exon][2]):
            if include_interchromosomal or (not read.is_paired or (read.is_paired and read.next_reference_name in list(set([_[0] for _ in exons.values()])))):
                if include_duplicates or (not read.is_duplicate):
                    for sj in get_splice_junction_positions(read):
                        if sj[0] != None and sj[1] == (exons['24'][1]):
                            
                            if sj[0] == exons['25'][2]:
                                set_sv1.add(read.query_name)
                                break
                            elif sj[0] == exons['26'][2]:
                                set_wt.add(read.query_name)
                                break


    exon = "25"
    for read in fh.fetch(exons[exon][0], exons[exon][1], exons[exon][2]):
        if read.get_overlap(exons[exon][1], exons[exon][2]):
            if include_interchromosomal or (not read.is_paired or (read.is_paired and read.next_reference_name in list(set([_[0] for _ in exons.values()])))):
                if include_duplicates or (not read.is_duplicate):
                    for sj in get_splice_junction_positions(read):
                        if sj[0] != None and sj[1] == (exons['25'][1]) and sj[0] == exons['26'][2]:
                            set_sv2.add(read.query_name)
                            break


    wt = set_wt.difference(set_sv1, set_sv2)
    
    splice_right_intron = set_sv1.difference(set_sv2, set_wt)
    splice_left_intron = set_sv2.difference(set_sv1, set_wt)
    splice_both_introns = set_sv1.intersection(set_sv2).difference(set_wt)
    
    discrepancies = set_wt.intersection(set_sv1.union(set_sv2))

    if discrepancies:
        print( "Warning, found (n="+str(len(discrepancies))+") read with SJ's to both wt as splice variant: " + list(discrepancies)[0], file=sys.stderr)

    
    return {'7B89: wt read [24 -> 26]': wt,
            '7B89: sv read [24 -> 25 -> 26]': splice_both_introns,
            '7B89: sv read [24 -> 25]': splice_right_intron,
            '7B89: sv read [25 -> 26]': splice_left_intron,
            '7B89: discrepant read': discrepancies} # this does not test for discrepancies



def extract_FN1_11A12_sv_reads_based_on_sjs(bam, exons, include_interchromosomal, include_duplicates):
    set_wt = set()# readnames of those that splice from exon 31 to 33
    set_sv1 = set()# readnames of those that splice from exon 31 to 32
    set_sv2 = set()# readnames of those that splice from exon 32 to 33
    
    fh = pysam.AlignmentFile(bam, "rb")
    exons = check_or_update_chr(exons, fh)

    exon = "31"
    for read in fh.fetch(exons[exon][0], exons[exon][1], exons[exon][2]):
        if read.get_overlap(exons[exon][1], exons[exon][2]):
            if include_interchromosomal or (not read.is_paired or (read.is_paired and read.next_reference_name in list(set([_[0] for _ in exons.values()])))):
                if include_duplicates or (not read.is_duplicate):
                    for sj in get_splice_junction_positions(read):
                        if sj[0] != None and sj[1] == (exons['31'][1]):
                            
                            if sj[0] == exons['32'][2]:
                                set_sv1.add(read.query_name)
                                break
                            elif sj[0] == exons['33'][2]:
                                set_wt.add(read.query_name)
                                break


    exon = "32"
    for read in fh.fetch(exons[exon][0], exons[exon][1], exons[exon][2]):
        if read.get_overlap(exons[exon][1], exons[exon][2]):
            if include_interchromosomal or (not read.is_paired or (read.is_paired and read.next_reference_name in list(set([_[0] for _ in exons.values()])))):
                if include_duplicates or (not read.is_duplicate):
                    for sj in get_splice_junction_positions(read):
                        if sj[0] != None and sj[1] == (exons['32'][1]) and sj[0] == exons['33'][2]:
                            set_sv2.add(read.query_name)
                            break


    wt = set_wt.difference(set_sv1, set_sv2)
    
    splice_right_intron = set_sv1.difference(set_sv2, set_wt)
    splice_left_intron = set_sv2.difference(set_sv1, set_wt)
    splice_both_introns = set_sv1.intersection(set_sv2).difference(set_wt)
    
    discrepancies = set_wt.intersection(set_sv1.union(set_sv2))

    if discrepancies:
        print( "Warning, found (n="+str(len(discrepancies))+") read with SJ's to both wt as splice variant: " + list(discrepancies)[0], file=sys.stderr)

    
    return {'11A12: wt read [31 -> 33]': wt,
            '11A12: sv read [31 -> 32 -> 33]': splice_both_introns,
            '11A12: sv read [31 -> 32]': splice_right_intron,
            '11A12: sv read [32 -> 33]': splice_left_intron,
            '11A12: discrepant read': discrepancies} # this does not test for discrepancies




