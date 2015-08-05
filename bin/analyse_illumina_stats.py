#!/usr/bin/env python
#
#     analyse_illumina_stats.py: summarise per lane/per sample stats
#     Copyright (C) University of Manchester 2015 Peter Briggs
#
#########################################################################
#
# analyse_illumina_stats.py
#
#########################################################################

"""analyse_illumina_stats

Utility to analyse the 'statistics.info' file from auto_process.py
and summarise the number of reads per sample in each lane (and the
percentage of reads that represents in the lane).

Example output:

Lane 1
Total reads = 182851745
- KatharineDibb/KD-K1	79888058	43.7%
- KatharineDibb/KD-K3	97854292	53.5%
- Undetermined_indices/lane1	5109395	2.8%
...

"""

#######################################################################
# Module metadata
#######################################################################

__version__ = "0.0.1"

#######################################################################
# Imports
#######################################################################

import optparse
import bcftbx.TabFile as tf
from bcftbx.IlluminaData import IlluminaFastq

#######################################################################
# Main script
#######################################################################

if __name__ == '__main__':

    # Process command line
    p = optparse.OptionParser(usage="%prog [OPTIONS] STATS_INFO_FILE",
                              version="%prog "+__version__,
                              description="Summarise the per sample stats for each lane "
                              "for an Illumina sequencing run, using the data from "
                              "STATS_INFO_FILE (typically called 'statistics.info' and "
                              "found in the top-level directory of sequencing runs "
                              "processed using auto_process.py.")
    options,args = p.parse_args()
    if len(args) != 1:
        p.error("expects a single argument (statistics file)")

    # Read in data
    stats = tf.TabFile(args[0],first_line_is_header=True)
    data = dict()
    for line in stats:
        # Collect sample name, lane etc
        fq = IlluminaFastq(line['Fastq'])
        if fq.read_number != 1:
            # Only interested in R1 reads
            continue
        lane = fq.lane_number
        sample = "%s/%s" % (line['Project'],line['Sample'])
        nreads = line['Nreads']
        # Update information in dict
        if lane not in data:
            data[lane] = dict()
        try:
            data[lane][sample] += nreads
        except KeyError:
            data[lane][sample] = nreads
    # Get list of lanes
    lanes = [int(x) for x in data.keys()]
    lanes.sort()
    # Report
    for lane in lanes:
        print "\nLane %d" % lane
        samples = data[lane].keys()
        samples.sort()
        total_reads = sum([data[lane][x] for x in samples])
        print "Total reads = %d" % total_reads
        for sample in samples:
            nreads = data[lane][sample]
            print "- %s\t%d\t%.1f%%" % (sample,nreads,float(nreads)/total_reads*100.0)
    
            
