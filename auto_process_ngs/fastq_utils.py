#!/usr/bin/env python
#
#     fastq_utils.py: utility functions for operating on fastq files
#     Copyright (C) University of Manchester 2016 Peter Briggs
#
########################################################################
#
# fastq_utils.py
#
#########################################################################

"""
fastq_utils.py

Utility functions for operating on Fastq files:

- collect_fastqs: return Fastqs from CASAVA/bcl2fastq output dir
- assign_barcodes_single_end: extract and assign inline barcodes

"""

#######################################################################
# Imports
#######################################################################

import os
import gzip
from bcftbx.IlluminaData import IlluminaData
from bcftbx.IlluminaData import IlluminaFastq
from bcftbx.IlluminaData import IlluminaDataError
from bcftbx.FASTQFile import FastqIterator

#######################################################################
# Functions
#######################################################################

def collect_fastqs(dirn,lane=None,read=None,unaligned_dir=None):
    """
    Automatically collect Fastq files from bcl2fastq/CASAVA output dir

    Arguments:
      dirn (str): path to the directory to collect files from (can
        be either the bcl2fastq/CASAVA output dir, or the directory
        above in which case the `unaligned_dir` should also be
        specified)
      lane (integer): optionally restrict the Fastqs to those with
        the specified lane in their names (default is to return Fastqs
        for all lanes)
      read_number (integer): optionally restrict the Fastqs to the
        specified read number (1 or 2) (default is to return both R1
        and R2)
      unaligned_dir (str): specify the subdirectory of `dirn` with
        the bcl2fastq/CASAVA outputs (if `dirn` is not that directory)

    Returns:
      List: either a list of R1/R2 Fastq filename pairs, or a "flat"
        list of either R1 or R2 Fastq filenames (if ``read`` argument
        was specified)

    """
    # Load data into IlluminaData instance
    if unaligned_dir is None:
        unaligned_dir = os.path.basename(dirn.rstrip(os.sep))
        dirn = os.path.dirname(os.path.abspath(dirn.rstrip(os.sep)))
    try:
        illumina_data = IlluminaData(dirn,unaligned_dir=unaligned_dir)
    except Exception as ex:
        raise IlluminaDataError("Unable to read fastqs from %s: %s\n" %
                                (dirn,ex))
    # Collect fastq R1/R2 pairs
    fastqs_r = dict()
    fastqs_r[1] = []
    fastqs_r[2] = []
    for project in illumina_data.projects:
        for sample in project.samples:
            fastqs_r[1].extend(sample.fastq_subset(read_number=1,
                                                   full_path=True))
            fastqs_r[2].extend(sample.fastq_subset(read_number=2,
                                                   full_path=True))
    if illumina_data.undetermined:
        for sample in illumina_data.undetermined.samples:
            fastqs_r[1].extend(sample.fastq_subset(read_number=1,
                                                   full_path=True))
            fastqs_r[2].extend(sample.fastq_subset(read_number=2,
                                                   full_path=True))
    # Filter by lane
    if lane is not None:
        fastqs_r[1] = filter(lambda fq: IlluminaFastq(fq).lane_number == lane,
                             fastqs_r[1])
        fastqs_r[2] = filter(lambda fq: IlluminaFastq(fq).lane_number == lane,
                             fastqs_r[2])
    # Sort into order
    fastqs_r[1].sort()
    fastqs_r[2].sort()
    # Handle when only one read is selected
    if read is not None:
        fastqs = fastqs_r[read]
    else:
        if not illumina_data.paired_end:
            fastqs = fastqs_r[1]
        else:
            fastqs = []
            for fq1,fq2 in zip(fastqs_r[1],fastqs_r[2]):
                fastqs.extend([fq1,fq2])
    return fastqs

def assign_barcodes_single_end(fastq_in,fastq_out,n=5):
    """
    Extract inline barcodes and assign to Fastq read headers

    Strips the first n bases from each read of the input
    FASTQ file and assigns it to the index sequence for that
    read in the output file.

    If the supplied output file name ends with '.gz' then it
    will be gzipped.

    Arguments:
      fastq_in (str): input FASTQ file (can be gzipped)
      fastq_out (str): output FASTQ file (will be gzipped if
        ending with '.gz')
      n (integer): number of bases to extract and assign as
        index sequence (default: 5)

    Returns:
      Integer: number of reads processed.

    """
    if fastq_out.endswith('.gz'):
        fp = gzip.GzipFile(filename=fastq_out,mode='wb')
    else:
        fp = open(fastq_out,'w')
    print "Processing reads from %s" % fastq_in
    nread = 0
    for read in FastqIterator(fastq_in):
        # Extract new barcode sequence
        barcode = read.sequence[:n]
        # Truncate sequence and quality accordingly
        sequence = read.sequence[n:]
        quality = read.quality[n:]
        # Assign new values and write to output
        read.seqid.index_sequence = barcode
        read.sequence = sequence
        read.quality = quality
        fp.write("%s\n" % read)
        nread += 1
    print "Finished (%d reads processed)" % nread
    return nread
