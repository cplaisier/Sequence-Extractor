#################################################################
# @Program: sequenceExtractor.py                                #
# @Version: 2                                                   #
# @Author: Chris Plaisier                                       #
# @Sponsored by:                                                #
# Nitin Baliga, ISB                                             #
# Institute for Systems Biology                                 #
# 1441 North 34th Street                                        #
# Seattle, Washington  98103-8904                               #
# (216) 732-2139                                                #
# @Also Sponsored by:                                           #
# Luxembourg Systems Biology Grant                              #
#                                                               #
# If this program is used in your analysis please mention who   #
# built it. Thanks. :-)                                         #
#                                                               #
# Copyrighted by Chris Plaisier  7/14/2010                      #
#################################################################

import os, gzip
from copy import deepcopy
from ftplib import FTP
from subprocess import *
import tarfile
import numpy

promoterSeq = [500,-200]
#min3pUTR = 831 # PMID = 11465035

# Convert exons into introns
def exon2intron(geneCoords):
    introns = []
    exonStarts = [int(x) for x in geneCoords['exonStarts'].split(',') if x]
    exonEnds = [int(x) for x in geneCoords['exonEnds'].split(',') if x]    
    exonStarts.pop(0)
    exonEnds.pop(-1)
    for intron in range(len(exonEnds)):
        introns.append([exonEnds[intron],exonStarts[intron]])
    return introns
    
# Complement
def complement(seq):
    complement = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N'}
    complseq = [complement[base] for base in seq]
    return complseq

# Reverse complement
def reverseComplement(seq):
    seq = list(seq)
    seq.reverse()
    return ''.join(complement(seq))

# Function to retreive boundaries for 3pUTR
def get3pUTR(geneCoords,min3pUTR):
    tmpUTR = []
    # Set boundaries
    if geneCoords['strand']=='+':
        start = geneCoords['cdsEnd']
        end = geneCoords['txEnd']
        tmpUTR = [[start,end]]
    elif geneCoords['strand']=='-':
        start = geneCoords['txStart']
        end = geneCoords['cdsStart']
        tmpUTR = [[start,end]]
    # Screen to see if introns exist in this region
    introns = exon2intron(geneCoords)
    for intron in range(len(introns)):
        for part in range(len(tmpUTR)):
            # If the exon lies in the 3' UTR region then take it out
            if tmpUTR[part][0] < introns[intron][0] < introns[intron][1] < tmpUTR[part][1]:
                tmpPart = tmpUTR.pop(part)
                # Insert in reverse order to preserve order of array
                tmpUTR.insert(part,[introns[intron][1],tmpPart[1]])
                tmpUTR.insert(part,[tmpPart[0],introns[intron][0]])
                # Now break from the loop back to the exon level loop to goto next exon
                break
    if lenICS(tmpUTR)<min3pUTR:
        diff = min3pUTR - lenICS(tmpUTR)
        if geneCoords['strand']=='+':
            tmpUTR[len(tmpUTR)-1][1] += diff
        elif geneCoords['strand']=='-':
            tmpUTR[0][0] = tmpUTR[0][0] - diff
    return tmpUTR

# Function to retreive boundaries for 5pUTR
def get5pUTR(geneCoords,min5pUTR):
    tmpUTR = []
    # Set boundaries
    if geneCoords['strand']=='+':
        start = geneCoords['txStart']
        end = geneCoords['cdsStart']
        tmpUTR = [[start,end]]
    elif geneCoords['strand']=='-':
        start = geneCoords['txEnd']
        end = geneCoords['cdsEnd']
        tmpUTR = [[start,end]]
    # Screen to see if introns exist in this region
    introns = exon2intron(geneCoords)
    for intron in range(len(introns)):
        for part in range(len(tmpUTR)):
            # If the exon lies in the 5' UTR region then take it out
            if tmpUTR[part][0] < introns[intron][0] < introns[intron][1] < tmpUTR[part][1]:
                tmpPart = tmpUTR.pop(part)
                # Insert in reverse order to preserve order of array
                tmpUTR.insert(part,[introns[intron][1],tmpPart[1]])
                tmpUTR.insert(part,[tmpPart[0],introns[intron][0]])
                # Now break from the loop back to the exon level loop to goto next exon
                break
    if lenICS(tmpUTR)<min5pUTR:
        diff = min5pUTR - lenICS(tmpUTR)
        if geneCoords['strand']=='+':
            tmpUTR[len(tmpUTR)-1][1] += diff
        elif geneCoords['strand']=='-':
            tmpUTR[0][0] = tmpUTR[0][0] - diff
    return tmpUTR

# Function to retreive boundaries for CDS
def getCDS(geneCoords):
    tmpCDS = []
    # Set boundaries
    start = geneCoords['cdsStart']
    end = geneCoords['cdsEnd']
    tmpCDS = [[start,end]]
    # Screen to see if introns exist in this region
    introns = exon2intron(geneCoords)
    for intron in range(len(introns)):
        for part in range(len(tmpCDS)):
            # If the exon lies in the 3' UTR region then take it out
            if tmpCDS[part][0] < introns[intron][0] < introns[intron][1] < tmpCDS[part][1]:
                tmpPart = tmpCDS.pop(part)
                # Insert in reverse order to preserve order of array
                tmpCDS.insert(part,[introns[intron][1],tmpPart[1]])
                tmpCDS.insert(part,[tmpPart[0],introns[intron][0]])
                # Now break from the loop back to the exon level loop to goto next exon
                break
    return tmpCDS

# Function to retreive boundaries for 3pUTR
def getPromoter(geneCoords,upstream):
    if geneCoords['strand']=='+':
        return [(geneCoords['txStart'] - upstream[0]), (geneCoords['txStart'] - upstream[1])]
    elif geneCoords['strand']=='-':
        return [(geneCoords['txEnd'] + upstream[1]), (geneCoords['txEnd'] + upstream[0])]

# Function to get the length of a 3pUTR
def lenICS(utrCoords):
    utrLen = 0
    for part in utrCoords:
        utrLen = utrLen + (part[1] - part[0])
    return utrLen

# Extract the unique elements
def uniquify(str1):
    tmp = []
    splitUp = str1.strip().split(',')
    for i in splitUp:
        if not i in tmp:
            tmp.append(i)
    return tmp

# Merge overlaps and give back sequences [[5pStart,5pEnd], [[3pStart,3pEnd], ... ]]
# !!! - Assumes that the mergeDem entries come from the same chromosome
def mergeSeqs(mergeDem,upstream,min5pUTR,min3pUTR):
    orig = mergeDem[0]
    strand = orig['strand']
    # Grab the starting promoter, 5' UTR, coding and 3' UTR sequences
    promoter = getPromoter(orig,upstream)
    p5utr = get5pUTR(orig,min5pUTR)
    cds = getCDS(orig)
    p3utr = get3pUTR(orig,min3pUTR)
    # Now iterate through the rest and merge
    for i in range(1,len(mergeDem)):
        mergeMe = mergeDem[i]
        # Merge promoter sequences
        promoterM = getPromoter(mergeMe,upstream)
        if promoterM[0] < promoter[0] <= promoterM[1] <= promoter[1]:
            promoter[0] = promoterM[0]
        elif promoter[0] <= promoterM[0] <= promoter[1] < promoterM[1]:
            promoter[1] = promoterM[1]
        elif not promoter==promoterM:
            # COMPROMISE HERE: Then will take the one for the longest transcript
            if strand=='-':
                if promoterM[0] > promoter[1]:
                    promoter = promoterM
            elif strand=='+':
                if promoterM[1] < promoter[0]:
                    promoter = promoterM
        # Merge 5' UTR
        p5utrM = get5pUTR(mergeMe,min5pUTR)
        #print p5utrM,'; ',p5utr
        if not ((p5utrM[0][0]==p5utr[0][0]) and (p5utrM[len(p5utrM)-1][1]==p5utr[len(p5utr)-1][1])):
            #print "Need merging!"
            #print p5utrM,p5utr
            p5utr = mergeICS(p5utr,p5utrM)
            #print "Merged: ",p3utr
        # Merge coding sequences
        cdsM = getCDS(mergeMe)
        #print cdsM,'; ',cds
        if not ((cdsM[0][0]==cds[0][0]) and (cdsM[len(cdsM)-1][1]==cds[len(cds)-1][1])):
            #print "Need merging!"
            #print cdsM,cds
            cds = mergeICS(cds,cdsM)
            #print "Merged: ",cds
        # Merge 3' UTR
        p3utrM = get3pUTR(mergeMe,min3pUTR)
        #print p3utrM,'; ',p3utr
        if not ((p3utrM[0][0]==p3utr[0][0]) and (p3utrM[len(p3utrM)-1][1]==p3utr[len(p3utr)-1][1])):
            #print "Need merging!"
            #print p3utrM,p3utr
            p3utr = mergeICS(p3utr,p3utrM)
            #print "Merged: ",p3utr
    return [promoter, p5utr, cds, p3utr]

# Merge intron containing sequence (5' UTR, CDS and 3'UTR)
def mergeICS(p3utr,p3utrM):
    # If the start of the new 3' UTR is further down
    if p3utrM[0][0] < p3utr[0][0]:
        # If the other bound is equal the current 3' UTR
        # then just add the new upper bound to the UTR 
        if (p3utrM[0][1] == p3utr[0][1]) or (p3utr[0][0] <= p3utrM[0][1] <= p3utr[0][1]):
            p3utr[0][0] = p3utrM[0][0]
        # If the first bounds don't overlap then maybe a new intron?
        elif p3utrM[0][1] < p3utr[0][0]:
            insertMe = [p3utrM[0]]
            # Now check out the rest of the segments till one overlaps
            # with the first of the original set
            if len(p3utrM)>1:
                for part in range(1,len(p3utrM)):
                    if p3utrM[part][1] < p3utr[0][0]:
                        insertMe.append(p3utrM[part])
                    elif (p3utrM[part][0] < p3utr[0][0]) and (p3utr[0][0] <= p3utrM[part][1] <= p3utr[0][1]):
                        p3utr[0][0] = p3utrM[part][0]
                        break
                    else:
                        break
            # Now attach the new elements to the front
            p3utr = insertMe + p3utr
    # If the end of the new 3' UTR is farther forward
    if p3utrM[len(p3utrM)-1][1] > p3utr[len(p3utr)-1][1]:
        # If the other bound is equal the current 3' UTR
        # then just add the new upper bound to the UTR 
        if p3utrM[len(p3utrM)-1][0] == p3utr[len(p3utr)-1][0] or (p3utr[len(p3utr)-1][0] <= p3utrM[len(p3utrM)-1][0] <= p3utr[len(p3utr)-1][1]):
            p3utr[len(p3utr)-1][1] = p3utrM[len(p3utrM)-1][1]
        # If the first bound don't overlap then maybe a new intron?
        elif p3utrM[len(p3utrM)-1][0] > p3utr[len(p3utr)-1][1]:
            appendMe = [p3utrM[0]]
            # Now check out the rest of the segments till one overlaps
            # with the first of the original set
            if len(p3utrM)>1:
                for part in range(1,len(p3utrM)):
                    if p3utrM[part][0] > p3utr[len(p3utr)-1][1]:
                        appendMe.insert(0,p3utrM[part])
                    elif (p3utrM[part][1] > p3utr[len(p3utr)-1][1]) and (p3utr[0][0] <= p3utrM[part][1] <= p3utr[0][1]):
                        p3utr[len(p3utr)-1][1] = p3utrM[part][1]
                        break
                    else:
                        break
            # Now append the new elements to the end
            p3utr = p3utr + appendMe
    return p3utr 

# Download gene identifier conversion table from NCBI if not already done
print 'Downloading converstion table for Entrez IDs to RefSeq IDs...'
if not os.path.exists('gene2refseq.gz'):
    ftp1 = FTP('ftp.ncbi.nih.gov')
    ftp1.login()
    ftp1.cwd('/gene/DATA/')
    outFile = open('gene2refseq.gz','wb')
    ftp1.retrbinary('RETR gene2refseq.gz',outFile.write)
    outFile.close()
    ftp1.quit()



# Start cycling through species need to have the orgId which can be had from NCBI taxonomy website, and the name of the sequence file from UCSC FTP site
# Right now only setup to use where the sequence files is a chromFaMasked.tar.gz in the bigZip directory. Otherwise can't do that species yet.
# Should be easy to setup up just not ready to screw around with that yet. Probably just read in from chromosomes dir and then convert lowercase to 'N'.
organismCodes = {}
# Get organisms three letter codes from miRBase
ftp1 = FTP('mirbase.org')
ftp1.login()
ftp1.cwd('/pub/mirbase/CURRENT')
outFile = open('organisms.txt','wb')
ftp1.retrbinary('RETR organisms.txt',outFile.write)
outFile.close()
ftp1.quit()
inFile = open('organisms.txt','r')
for line in inFile.readlines():
    splitUp = line.strip().split('\t')
    organismCodes[splitUp[2]] = { 'miRBase_3lc':splitUp[0] }

inFile.close()

# Get taxonomic ids for all species from NCBI
ftp1 = FTP('ftp.ncbi.nih.gov')
ftp1.login()
ftp1.cwd('/pub/taxonomy')
outFile = open('taxdump.tar.gz','wb')
ftp1.retrbinary('RETR taxdump.tar.gz',outFile.write)
outFile.close()
ftp1.quit()
tar = tarfile.open('taxdump.tar.gz')
tar.extractall()
tar.close()
inFile = open('names.dmp','r')
cnt = 0
for line in inFile.readlines():
    splitUp = line.strip().split('\t')
    if splitUp[2] in organismCodes:
        organismCodes[splitUp[2]]['NCBI_taxonId'] = splitUp[0]
        cnt += 1

inFile.close()
#print cnt

ftp1 = FTP('hgdownload.cse.ucsc.edu')
ftp1.login()
ftp1.cwd('/goldenPath/currentGenomes')
species = ftp1.nlst()
species.remove('.')
species.remove('..')
orgDict = {}
for organism in species:
    org1 = organism.replace('_',' ')
    if org1 in organismCodes:
        ftp1.cwd('/goldenPath/currentGenomes/'+organism+'/bigZips/')
        files1 = ftp1.nlst()
        ftp1.cwd('/goldenPath/currentGenomes/'+organism+'/database/')
        files2 = ftp1.nlst()
        if 'chromFaMasked.tar.gz' in files1 and 'refGene.txt.gz' in files2:
            print org1
            orgDict[organism] = { 'orgId':int(organismCodes[org1]['NCBI_taxonId']), 'seqFile': 'chromFaMasked.tar.gz' }

#orgDict = { 'Drosophila_melanogaster': {'orgId':7227, 'seqFile':'chromFaMasked.tar.gz'} }
#orgDict = { 'Homo_sapiens': {'orgId':9606, 'seqFile':'chromFaMasked.tar.gz'}, 'Drosophila_melanogaster': {'orgId':7227, 'seqFile':'chromFaMasked.tar.gz'}, 'Gallus_gallus': {'orgId':9031, 'seqFile':'chromFaMasked.tar.gz'}, 'Mus_musculus': {'orgId':10090, 'seqFile':'chromFaMasked.tar.gz'} } #, 'Caenorhabditis_elegans':{'orgId':6239, 'seqFile':'chromFaMasked.tar.gz'}, 'Pan_troglodytes':{'orgId':9598, 'seqFile':'chromFaMasked.tar.gz'}, 'Rattus_norvegicus':{'orgId':10116, 'seqFile':'chromFaMasked.tar.gz'}, 'Pongo_pygmaeus_abelii':{'orgId':9601, 'seqFile':'chromFaMasked.tar.gz'}, 'Canis_familiaris':{'orgId':9615, 'seqFile':'chromFaMasked.tar.gz'}, 'Rhesus_macaque':{'orgId':9544, 'seqFile':'chromFaMasked.tar.gz'}, 'Equus_caballus':{'orgId':9796, 'seqFile':'chromFaMasked.tar.gz'}, 'Monodelphis_domestica':{'orgId':13616, 'seqFile':'chromFaMasked.tar.gz'}, 'Taeniopygia_guttata':{'orgId':59729, 'seqFile':'chromFaMasked.tar.gz'} } #Can add others upon request
#orgDict = { 'Caenorhabditis_elegans': {'orgId':6239, 'seqFile':'chromFaMasked.tar.gz'} }
#orgDict = { 'Canis_familiaris':{'orgId':9615, 'seqFile':'chromFaMasked.tar.gz'}, 'Rattus_norvegicus':{'orgId':10116, 'seqFile':'chromFaMasked.tar.gz'} }
#orgDict = { 'Rattus_norvegicus':{'orgId':10116, 'seqFile':'chromFaMasked.tar.gz'} }
#orgDict = { 'Homo_sapiens': {'orgId':9606, 'seqFile':'chromFaMasked.tar.gz'} }
#orgDict = { 'Homo_sapiens': {'orgId':9606, 'seqFile':'chromFaMasked.tar.gz'}, 'Drosophila_melanogaster': {'orgId':7227, 'seqFile':'chromFaMasked.tar.gz'}, 'Gallus_gallus': {'orgId':9031, 'seqFile':'chromFaMasked.tar.gz'}, 'Mus_musculus': {'orgId':10090, 'seqFile':'chromFaMasked.tar.gz'}, 'Caenorhabditis_elegans':{'orgId':6239, 'seqFile':'chromFaMasked.tar.gz'}, 'Pan_troglodytes':{'orgId':9598, 'seqFile':'chromFaMasked.tar.gz'}, 'Rattus_norvegicus':{'orgId':10116, 'seqFile':'chromFaMasked.tar.gz'}, 'Pongo_pygmaeus_abelii':{'orgId':9601, 'seqFile':'chromFaMasked.tar.gz'}, 'Canis_familiaris':{'orgId':9615, 'seqFile':'chromFaMasked.tar.gz'}, 'Rhesus_macaque':{'orgId':9544, 'seqFile':'chromFaMasked.tar.gz'}, 'Equus_caballus':{'orgId':9796, 'seqFile':'chromFaMasked.tar.gz'}, 'Monodelphis_domestica':{'orgId':13616, 'seqFile':'chromFaMasked.tar.gz'}, 'Taeniopygia_guttata':{'orgId':59729, 'seqFile':'chromFaMasked.tar.gz'} } #Can add others upon request
#orgDict = { 'Homo_sapiens': {'orgId':9606, 'seqFile':'chromFaMasked.tar.gz'}, 'Drosophila_melanogaster': {'orgId':7227, 'seqFile':'chromFaMasked.tar.gz'}, 'Gallus gallus': {'orgId':9031, 'seqFile':'chromFaMasked.tar.gz'}, 'Mus_musculus': {'orgId':10090, 'seqFile':'chromFaMasked.tar.gz'}, 'Danio_rerio': {'orgId':7955, 'seqFile':'danRer6.fa.masked.gz'} } # This may be useful when I figure out how to get the sequence from zebrafish and cow
#orgDict = { 'Homo_sapiens': {'orgId':9606, 'seqFile':'chromFaMasked.tar.gz'} }
# 'Sus_scrofa':{'orgId':9823, 'seqFile':'chromFaMasked.tar.gz'} - pig doesn't have refGene.txt.gz
# 'Takifugu_rubripes':{'orgId':31033, 'seqFile':'chromFaMasked.tar.gz'} - doesn't have refGene.txt.gz
# 'Tetraodon_nigroviridis':{'orgId':99883, 'seqFile':'chromFaMasked.tar.gz'} - doesn't have refGene.txt.gz
orgData = {}
for org in orgDict:
    print 'Starting on '+str(org)+'...'
    print '  Downloading genomic data...'
    if not os.path.exists(str(org)+'/chrs'):
        os.makedirs(str(org)+'/chrs')
    # Download genome information for organism if not already done
    if not os.path.exists(str(org)+'/chrs/'+str(orgDict[org]['seqFile'])) or not os.path.exists(str(org)+'/refGene.txt.gz'):    
        # Download genome information from UCSC
        ftp1 = FTP('hgdownload.cse.ucsc.edu')
        ftp1.login()
        # Get the chromosome data
        ftp1.cwd('/goldenPath/currentGenomes/'+str(org)+'/bigZips/')
        outFile = open(str(org)+'/chrs/'+str(orgDict[org]['seqFile']),'wb')
        ftp1.retrbinary('RETR '+str(orgDict[org]['seqFile']),outFile.write)
        outFile.close()
        # Get gene data
        ftp1.cwd('/goldenPath/currentGenomes/'+str(org)+'/database/')
        outFile = open(str(org)+'/refGene.txt.gz','wb')
        ftp1.retrbinary('RETR refGene.txt.gz',outFile.write)
        outFile.close()
        ftp1.quit()

    print '  Building dictionaries...'
    # 1. Read in refSeqCoords
    inFile = gzip.open(str(org)+'/refGene.txt.gz','r')
    refseqCoords = {}
    chrs = []
    while 1:
        line = inFile.readline()
        if not line:
            break
        splitUp = line.strip().split('\t')
        if len(splitUp)>=13 and splitUp[13]=='cmpl':
            if not splitUp[1] in refseqCoords:
                if not splitUp[2] in chrs:
                    chrs.append(splitUp[2])
                refseqCoords[splitUp[1]] = {'chr':splitUp[2], 'strand':splitUp[3], 'txStart':int(splitUp[4]), 'txEnd':int(splitUp[5]), 'cdsStart':int(splitUp[6]), 'cdsEnd':int(splitUp[7]), 'exonCount':int(splitUp[8]), 'exonStarts':splitUp[9], 'exonEnds':splitUp[10], 'geneName':splitUp[12], 'exonFrames':[int(x) for x in splitUp[15].split(',') if x]}
        elif not len(splitUp)>=13:
            if not splitUp[0] in refseqCoords:
                if not splitUp[1] in chrs:
                    chrs.append(splitUp[1])
                # Build the exonFrames determine which 
                refseqCoords[splitUp[0]] = {'chr':splitUp[1], 'strand':splitUp[2], 'txStart':int(splitUp[3]), 'txEnd':int(splitUp[4]), 'cdsStart':int(splitUp[5]), 'cdsEnd':int(splitUp[6]), 'exonCount':int(splitUp[7]), 'exonStarts':splitUp[8], 'exonEnds':splitUp[9]}
    inFile.close()

    # 1a. Calculate the median 5' & 3' UTR length to be used for unknowns
    p5utrLens = []
    p3utrLens = []
    for refseq in refseqCoords:
        a1 = lenICS(get5pUTR(refseqCoords[refseq],0))
        if not a1==0:
            p5utrLens.append(a1)
        a1 = lenICS(get3pUTR(refseqCoords[refseq],0))
        if not a1==0:
            p3utrLens.append(a1)
    min5pUTR = int(numpy.median(p5utrLens))
    min3pUTR = int(numpy.median(p3utrLens))
    orgData[org] = { 'min5pUTR': min5pUTR, 'min3pUTR':min3pUTR }
    print '  Median 5\' UTR length =',min5pUTR,'bp (n =',len(p5utrLens),')'
    print '  Median 3\' UTR length =',min3pUTR,'bp (n =',len(p3utrLens),')'

    # 2. Make a dictionary of EntrezIDs to RefSeqIds
    inFile = gzip.open('gene2refseq.gz','r')
    inFile.readline() # skip header
    entrezId2refSeq = {}
    while 1:
        line = inFile.readline()
        if not line:
            break
        # Only add those that have the correct NCBI organism ID
        splitUp = line.strip().split('\t')
        if int(splitUp[0])==orgDict[org]['orgId']:
            #print splitUp[3],splitUp[3].split('.')[0]
            # Check that the nucleotide ID is not a '-' and that it has genomic coordiantes assocaited with it
            if not splitUp[3]=='-' and splitUp[3].split('.')[0] in refseqCoords:
                if not int(splitUp[1]) in entrezId2refSeq:
                    entrezId2refSeq[int(splitUp[1])] = [splitUp[3].split('.')[0]]
                else:
                    entrezId2refSeq[int(splitUp[1])].append(splitUp[3].split('.')[0])
    inFile.close()
    print ' ',len(entrezId2refSeq),len(refseqCoords)
    orgData[org]['entrez2refSeq'] = len(entrezId2refSeq)
    orgData[org]['refseqCoords'] = len(refseqCoords)
    
    print '  Now collapsing and merging RefSeq IDs into Entrez IDs...'
    # 3. Merege multiple refseq IDs corresponding to a single entrezID
    # Now merge the data
    #chrs = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY'] # This needs to come from the chromsome fasta masked files
    mergedSet = {}
    for chr in chrs:
        mergedSet[chr] = {}
    baddies = []
    chrNoMatch = 0
    gotGenes = []
    for entrezId in entrezId2refSeq:
        chr = refseqCoords[entrezId2refSeq[entrezId][0]]['chr']
        strand = refseqCoords[entrezId2refSeq[entrezId][0]]['strand']
        if chr in mergedSet:
            if len(entrezId2refSeq[entrezId])>1:
                # There are duplicates so build a list of their refseqCoords and merge them
                mergeDem = []
                for refseq in entrezId2refSeq[entrezId]:
                    tot = 0
                    negOne = 0
                    goodOrBad = 1
                    if 'exonFrames' in refseqCoords[refseq]:
                        for i in refseqCoords[refseq]['exonFrames']:
                            if i == -1:
                                negOne += 1
                            tot += 1
                        if len(refseqCoords[refseq]['exonFrames'])>=5:
                            goodOrBad = 1-float(negOne)/float(tot)
                    #print goodOrBad
                    if chr=='' and goodOrBad>=0.5:
                        chr = refseqCoords[refseq]['chr']
                        mergeDem.append(refseqCoords[refseq])
                    elif chr==refseqCoords[refseq]['chr'] and refseqCoords[refseq]['strand'] and goodOrBad>=0.5:
                        mergeDem.append(refseqCoords[refseq])
                    elif goodOrBad<0.5:
                        #print 'Baddie taken out: ',refseq
                        baddies.append(refseq)
                    else:
                        #print 'Uh oh! Chr and Strand don\'t match! EntrezID = ',entrezId,'; refseqID = ',refseq
                        chrNoMatch += 1
                if len(mergeDem)>1:
                    mergedSet[chr][entrezId] = mergeSeqs(mergeDem,promoterSeq,min5pUTR,min3pUTR) + [strand]
                    gotGenes.append(entrezId)
                    #print entrez, len(entrez2refseq[entrez]), lenICS(mergedSet[chr][entrez][1])
                else:
                    promoter = getPromoter(refseqCoords[(entrezId2refSeq[entrezId])[0]],promoterSeq)
                    p5utr = get5pUTR(refseqCoords[(entrezId2refSeq[entrezId])[0]],min5pUTR)
                    cds = getCDS(refseqCoords[(entrezId2refSeq[entrezId])[0]])
                    p3utr = get3pUTR(refseqCoords[(entrezId2refSeq[entrezId])[0]],min3pUTR)
                    mergedSet[chr][entrezId] = [promoter,p5utr,cds,p3utr,strand]
                    gotGenes.append(entrezId)
            else:
                promoter = getPromoter(refseqCoords[(entrezId2refSeq[entrezId])[0]],promoterSeq)
                p5utr = get5pUTR(refseqCoords[(entrezId2refSeq[entrezId])[0]],min5pUTR)
                cds = getCDS(refseqCoords[(entrezId2refSeq[entrezId])[0]])
                p3utr = get3pUTR(refseqCoords[(entrezId2refSeq[entrezId])[0]],min3pUTR)
                mergedSet[chr][entrezId] = [promoter,p5utr,cds,p3utr,strand]
                gotGenes.append(entrezId)
    badFile = open(str(org)+'/baddies.txt','w')
    badFile.write('\n'.join(baddies))
    badFile.close()
    
    # 4. Get the sizes for the regulatory regions
    mergedSetFile = open(str(org)+'/mergedSetsLengthsRefSeq.csv','w')
    mergedSetFile.write('EntrezID,Promoter Length,5pUTR Length,CDS Length,3pUTR Length')
    for chr in mergedSet:
        for entrezId in mergedSet[chr]:
            promoterL = mergedSet[chr][entrezId][0][1] - mergedSet[chr][entrezId][0][0]
            p5utrL = lenICS(mergedSet[chr][entrezId][1])
            cdsL = lenICS(mergedSet[chr][entrezId][2])
            p3utrL = lenICS(mergedSet[chr][entrezId][3])
            mergedSetFile.write('\n'+str(entrezId)+','+str(promoterL)+','+str(p5utrL)+','+str(cdsL)+','+str(p3utrL))
    mergedSetFile.close() 
    
    print '  Extracting the sequence data...'
    # 5. Unzip sequences for extraction
    #errOut = open(str(org)+'/stderr.out','w')
    tar = tarfile.open(str(org)+'/chrs/'+str(orgDict[org]['seqFile']))
    tar.extractall(path=str(org)+'/chrs')
    tar.close()

    # 6. Extract the sequences
    promoterFile = open(str(org)+'/promoterSeqs_'+str(org)+'.csv','w')
    p5utrFile = open(str(org)+'/p5utrSeqs_'+str(org)+'.csv','w')
    cdsFile = open(str(org)+'/cdsSeqs_'+str(org)+'.csv','w')
    p3utrFile = open(str(org)+'/p3utrSeqs_'+str(org)+'.csv','w')
    for chr in mergedSet:
        if os.path.exists(str(org)+'/chrs/'+str(chr)+'.fa.masked'):
            chrSeqFile = open(str(org)+'/chrs/'+str(chr)+'.fa.masked','r')
        elif os.path.exists(str(org)+'/chrs/'+str(chr).lstrip('chr').replace('_random','')+'/'+str(chr)+'.fa.masked'):
                chrSeqFile = open(str(org)+'/chrs/'+str(chr).lstrip('chr').replace('_random','')+'/'+str(chr)+'.fa.masked','r')
        else:
            print 'FATAL ERROR!!!! Arghhh',chr,'(',str(chr).lstrip('chr').replace('_random',''),')does not have a seqeunce file!'
            break
        chrSeqFile.readline() # Get rid of header
        chrSeq = [x.strip().upper() for x in chrSeqFile.readlines()]
        chrSeq = ''.join(chrSeq)
        for entrezId in mergedSet[chr]:
            promSeq = chrSeq[(mergedSet[chr][entrezId][0][0]-1):(mergedSet[chr][entrezId][0][1]-1)]
            p5utrSeq = ''
            for part in mergedSet[chr][entrezId][1]:
                p5utrSeq += chrSeq[(part[0]-1):(part[1]-1)]
            cdsSeq = ''
            for part in mergedSet[chr][entrezId][2]:
                cdsSeq += chrSeq[(part[0]-1):(part[1]-1)]
            p3utrSeq = ''
            for part in mergedSet[chr][entrezId][3]:
                p3utrSeq += chrSeq[(part[0]-1):(part[1]-1)]
            if mergedSet[chr][entrezId][4] == '-':
                promSeq = reverseComplement(promSeq)
                p5utrSeq = reverseComplement(p5utrSeq)
                cdsSeq = reverseComplement(cdsSeq)
                p3utrSeq = reverseComplement(p3utrSeq)
            promoterFile.write(str(entrezId)+','+str(promSeq)+'\n')
            p5utrFile.write(str(entrezId)+','+str(p5utrSeq)+'\n')
            cdsFile.write(str(entrezId)+','+str(cdsSeq)+'\n')
            p3utrFile.write(str(entrezId)+','+str(p3utrSeq)+'\n')
    promoterFile.close()
    p5utrFile.close()
    cdsFile.close()
    p3utrFile.close()

    # 7. Cleanup the seqeunce data
    #files = os.listdir(str(org)+'/chrs/')
    #for file in files:
    #    if file.find('.fa.maksed'):
    #        os.remove(str(org)+'/chrs/'+str(file))
    #errOut = open(str(org)+'/stderr.out','w')
    #rmProc = Popen('rm '+str(org)+'/chrs/*.fa.masked')
    #output = weederProc.communicate()
    #errOut.close()
    print 'Done!\n'

#outFile = open('summary.csv','w')
#outFile.write('org,min5pUTR (bp),min3pUTR (bp)\n')
#outFile.write('\n'.join([','.join([org, organismCodes[organism.replace('_',' ')]['NCBI_taxonId'], organismCodes[organism.replace('_',' ')]['miRBase_3lc'], orgData[org]['min5pUTR'], orgData[org]['min3pUTR']]) for org in orgData]))
#outFile.close()

