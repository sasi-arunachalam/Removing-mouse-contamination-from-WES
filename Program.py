
#!/usr/bin/env python

import re
import os
import sys
import argparse
import pysam

def main():
    #Parse command line
    parser = argparse.ArgumentParser(description="This script selects mouse or human alignment")
    parser.add_argument("mouseAlignment",help="This file should be a alignment generated from novoalign")
    parser.add_argument("humanAlignment")
    parser.add_argument("outputPrefix")
    
    args = parser.parse_args()
    
    #containers
    mouseAlignments = {}
    humanAlignments = {}
    observedReads = {}
    
    #Read in mouse alignments
    print "Reading mouse alignments"
    samFile1 = pysam.Samfile(args.mouseAlignment,"rb")
    #outFileH = pysam.Samfile(args.outputPrefix + ".human.bam","wb",template=samFile1)
    counter = 0
    for read in samFile1:
        if not read.is_unmapped and not read.is_supplementary:
            counter += 1
            if counter % 1000000 == 0:
                print counter
            qual = read.mapping_quality
            if read.query_name not in mouseAlignments:
                mouseAlignments[read.query_name] = [-1,-1]
            idx = 1
            if read.is_read1:
                idx = 0
            a_score = None
            for t in read.get_tags():
                if t[0] == "AS":
                    a_score = t[1]
            mouseAlignments[read.query_name][idx] = a_score
            observedReads[read.query_name] = ""
    samFile1.close()
    
    #Read in mouse alignments
    print "Reading human alignments"
    samFile2 = pysam.Samfile(args.humanAlignment,"rb")
    #outFileM = pysam.Samfile(args.outputPrefix + ".mouse.bam","wb",template=samFile2)
    counter = 0
    for read in samFile2:
        if not read.is_unmapped:
            counter += 1
            if counter % 1000000 == 0:
                print counter
            qual = read.mapping_quality
            if read.query_name not in humanAlignments:
                humanAlignments[read.query_name] = [-1,-1]
            idx = 1
            if read.is_read1:
                idx = 0
            a_score = None
            for t in read.get_tags():
                if t[0] == "AS":
                    a_score = t[1]
            humanAlignments[read.query_name][idx] = a_score
            observedReads[read.query_name] = ""
    samFile2.close()
    
     
    #Counters
    humanOnly = 0
    mouseOnly = 0
    humanBetterBoth = 0
    mouseBetterBoth = 0
    humanBetterOne = 0
    mouseBetterOne = 0
    mixNoneBetter = 0
    mixHumanBetter = 0
    mixMouseBetter = 0
    tie = 0
    
    mouseReads = {}
    humanReads = {}
    tieReads = {}
    
#Create output files
        
    for read in observedReads.keys():
        if read in humanAlignments:
            if read in mouseAlignments: #Read in both, here is the hard part
                h1,h2 = humanAlignments[read]
                m1,m2 = mouseAlignments[read]
                if h1 > m1:
                    if h2 > m2:
                        humanBetterBoth += 1
                        humanReads[read] = "HB"
                    elif h2 == m2:
                        humanBetterOne += 1
                        humanReads[read] = "HO"
                    else:
                        code = None
                        if (h1-m1) > (m2-h2):
                            mixHumanBetter += 1
                            humanReads[read] = "HM"
                            code = "HM"
                        elif (h1-m1) < (m2-h2):
                            mixMouseBetter += 1
                            mouseReads[read] = "MH"
                            code = "MH"
                        else:
                            mixNoneBetter += 1
                elif h1 < m1:
                    if h2 < m2:
                        mouseBetterBoth += 1
                        mouseReads[read] = "MB"
                    elif h2 == m2:
                        mouseBetterOne += 1
                        mouseReads[read] = "MO"
                    else:
                        code = None
                        if (m1-h1) > (h2-m2):
                            mixMouseBetter += 1
                            mouseReads[read] = "MH"
                            code = "MH"
                        elif (m1-h1) < (h2-m2):
                            mixHumanBetter += 1
                            humanReads[read] = "HM"
                            code = "HM"
                        else:
                            mixNoneBetter += 1
                else:
                    if h2 > m2:
                        humanBetterOne += 1
                        humanReads[read] = "HO"
                    elif h2 < m2:
                        mouseBetterOne += 1
                        mouseReads[read] = "MO"
                    else:
                        humanReads[read] = "HT"
                        tie += 1
                        tieReads[read] = ""
            else:
                humanOnly += 1
                humanReads[read] = "H"
        elif read in mouseAlignments:
            mouseOnly += 1
            mouseReads[read] = "M"
        else:
            print "Should not be here, exiting"
            sys.exit(1)
    
    samFile1 = pysam.Samfile(args.humanAlignment,"rb")
    outFileH1 = pysam.Samfile(args.outputPrefix + "_human_all.bam","wb",template=samFile1)
    outFileH2 = pysam.Samfile(args.outputPrefix + "_human_H_HB.bam","wb",template=samFile1)
    outFileH3 = pysam.Samfile(args.outputPrefix + "_human_H_HB_HO.bam","wb",template=samFile1)
    for read in samFile1:
        if read.query_name in humanReads:
            tag = humanReads[read.query_name]
            read.set_tag("HC",tag)
            outFileH1.write(read)
            if tag in ["H","HO","HB"]:
                outFileH3.write(read)
                if tag in ["H","HB"]:
                    outFileH2.write(read)
            
            
            
    outFileH1.close()
    outFileH2.close()
    outFileH3.close()
    samFile1.close()
 samFile2 = pysam.Samfile(args.mouseAlignment,"rb")
    outFileM = pysam.Samfile(args.outputPrefix + "_mouse.bam","wb",template=samFile2)
    for read in samFile2:
        if read.query_name in mouseReads:
            read.set_tag("HC",mouseReads[read.query_name])
            outFileM.write(read)
    outFileM.close()
    samFile2.close()
    
    ofStats = open(args.outputPrefix + "_stats.txt","w")
    ofStats.write("Total_reads {0}\n".format(len(observedReads.keys())))
    ofStats.write("Human_only {0}\n".format(humanOnly))
    ofStats.write("Human_better {0}\n".format(humanBetterBoth))
    ofStats.write("Human_better_one {0}\n".format(humanBetterOne))
    ofStats.write("Human_better_mixed {0}\n".format(mixHumanBetter))
    ofStats.write("Mouse_only {0}\n".format(mouseOnly))
    ofStats.write("Mouse_better {0}\n".format(mouseBetterBoth))
    ofStats.write("Mouse_better_one {0}\n".format(mouseBetterOne))
    ofStats.write("Mouse_better_mixed {0}\n".format(mixMouseBetter))
    ofStats.write("Mixed, can't decide {0}\n".format(mixNoneBetter))
    ofStats.write("Tie {0}\n".format(tie))
    
    
if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print "user interrupted, exiting"
        sys.exit(1)

