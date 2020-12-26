#! /usr/bin/env python

from __future__ import print_function
from optparse import OptionParser
import sys
import pysam
import regex
import argparse
import os
import re
import collections
#from multiprocessing import Process, Queue, current_process, freeze_support
#import multiprocessing
import subprocess
import multiprocessing as mp
from multiprocessing import Lock
from multiprocessing import Manager

## global variable. get value from input
readseqlength = 0
ampliconfilecol = 9
## debug variable
debugPreProcessofSoftclip = 0
debugDetailedMode = 0
debugSimpleMode = 0

'''seqlen is the actual length of the read sequence, it is not included deletion, only include match, mismatch and insertion;
  reflen is the covered reference base according to reference genome;
  reflen is not included insertion, only include match, mismatch and deletion;
  cigar contains match/mismatch, insertion and deletion;
  MD contains match, mismatch, deletion;
  so, cigarlen - deletion = seqlen; cigarlen - insertion = MDlen;
'''

## MD='140T5A6^TT4' '9T130A5^T10'  -> ['140','T','5','A','6','^','T','T','4']
mdTagRe = re.compile('[0-9]+|\^+|[A-Z]')

## get DNA-seq`s reverse-complementary seq
def DNA_revcplm(sequence):
        return DNA_reverse(DNA_complement(sequence))

def DNA_reverse(sequence):
        sequence = sequence.upper()
        return sequence[::-1]

def DNA_complement(sequence):
        sequence = sequence.upper()
        sequence = sequence.replace('A', 't')
        sequence = sequence.replace('T', 'a')
        sequence = sequence.replace('C', 'g')
        sequence = sequence.replace('G', 'c')
        return sequence.upper()

def GetOutputFastqName(bamf):
        spltmp = filter(None, bamf.split('/'))
        fqpfx = spltmp[len(spltmp)-1].split('.')[0]
        return fqpfx

def GetBamFilePrefix(bamf):
        spltbf = filter(None, bamf.split('/'))
        bamfpfx = spltbf[len(spltbf)-1].split('.')[0]
        return bamfpfx

def MakeDir(path):
	#path=path.strip()
	#path=path.rstrip("\\")

	isExists=os.path.exists(path)
	if not isExists:
		os.makedirs(path)
		return True
	else:
		print('outpath already exists!')
		return False

def GetReadMapStrand(read):
        strd = ''
        if(read.is_reverse):
                strd = '-'
        else:
                strd = '+'
        return strd

def IsBamFile(bamname):
        bnl = bamname.split('.')
        if(bnl[len(bnl) - 1] == 'bam'):
                return 1
        else:
                return 0

class ReadSimpleInfo:
	def __init__(self):
		self.qname = ''
		self.seq = ''
		self.qual = ''
		self.chrom = ''
		self.mapspos = 0
		self.nextmapspos = 0
		self.strand = '+'
		self.cigar = [(0,readseqlength)]
		self.tlen = readseqlength
		self.NM = readseqlength
		self.MD = ''
		self.AS = readseqlength
		self.XS = 1
		self.mapepos = 0

	def update(self, read):
		self.qname = read.qname
		self.seq = read.query_sequence
		self.qual = read.query_qualities
		self.chrom = read.reference_name
		## in pysam, map start pos is 0-base position
		self.mapspos = read.reference_start + 1
		self.nextmapspos = read.next_reference_start + 1
		self.cigar = read.cigartuples
		self.tlen = read.template_length
		self.NM = read.get_tag('NM')
		self.MD = mdTagRe.findall(read.get_tag('MD'))
		self.AS = read.get_tag('AS')
		self.XS = read.get_tag('XS')
		## map end pos is the last base postion in the read  bug
		self.mapepos = read.reference_end
		
		if(not(read.is_reverse)):
			self.strand = '+'
		else:
			self.strand = '-'

	def modifyNM(self):
		## get del count from cigar
		delcnt = 0
		for x in self.cigar:
			## deletion
			if(x[0] == 2):
				delcnt += x[1]

		## get mismatch count from MD
		mismatchcnt = 0
		index = 0
		md = self.MD	
		while(index < len(md)):
			if(md[index].isalpha()):
				mismatchcnt += 1
			index += 1

		self.NM = delcnt + mismatchcnt

## seqlen: length of sequence in the read
class SeqInfoInRead(object):
	def __init__(self, mapstartpos, mapendpos, ampl_inss, ampl_inse):
		'''fp info
		'''
		if ampl_inss > mapstartpos:
			self.fpreflen = ampl_inss - mapstartpos
		else:
			self.fpreflen = 0
		self.fpseqlen = self.fpreflen
		self.fpinscnt = 0
		self.fpdelcnt = 0

		'''rp info
		'''
		if mapendpos > ampl_inse:
			self.rpreflen = mapendpos - ampl_inse
		else:
			self.rpreflen = 0
		self.rpseqlen = self.rpreflen
		self.rpinscnt = 0
		self.rpdelcnt = 0

	def updateFPSeqLen(self, cigar):
		if self.fpreflen == 0:
			return

		refbaselen = 0
		i = -1
		cg = cigar

		for x in cg:
			## cigar has M
			if(x[0] == 0):
				refbaselen += x[1]
			## cigar has I
			if(x[0] == 1):
				self.fpinscnt += x[1]
			## cigar has D, 90M2D30M
			if(x[0] == 2):
				refbaselen += x[1]
				self.fpdelcnt += x[1]
			i += 1

			if refbaselen >= self.fpreflen:
				break
		'''normal case: 103M2D50M

		'''
		if refbaselen == self.fpreflen:
			'''3 special case: fplen=20,   10M1I3M2D5M 1I90M       10M1I3M2D5M 1D90M        10M1I8M2D 5M1I90M
			'''
			if(cg[i + 1][0] == 1):
				self.fpinscnt += cg[i + 1][1]
			if(cg[i + 1][0] == 2):
				self.fpdelcnt += cg[i + 1][1]
				self.fpreflen += self.fpdelcnt

			if(debugDetailedMode == 1): 
				print("specialcase: fp refbaselen==self.fpreflen")
		else:
			'''bug:  1 special case.    fplen=20,    cigar:  15M8D 80M => reflen=13,seqlen=15
			'''
			if cg[i][0] == 2:
				self.fpreflen = refbaselen

				if(debugDetailedMode == 1): 
					print("specialcase: fp refbaselen > self.fpreflen and the count of deletion base is large")

		self.fpseqlen = self.fpreflen + self.fpinscnt - self.fpdelcnt

	def updateRPSeqLen(self, cigar):
		if self.rpreflen <= 0:
			return

		refbaselen = 0
		i = -1
		cg = cigar[::-1]

		for x in cg:
			## cigar has M
			if(x[0] == 0):
				refbaselen += x[1]
			## cigar has I
			if(x[0] == 1):
				self.rpinscnt += x[1]
			## cigar has D, 90M2D30M
			if(x[0] == 2):
				refbaselen += x[1]
				self.rpdelcnt += x[1]
			i += 1

			if(refbaselen >= self.rpreflen):
				break

		if refbaselen == self.rpreflen:
#			'''3 case: rplen=20, 90M1I 10M1I3M2D5M; fplen=20, 90M1D 10M1I3M2D5M; fplen=15,  90M1I 10M1I3M2D
#					 	i				 i 				    i	
#			'''
			if(cg[i + 1][0] == 1):
				self.rpinscnt += cg[i + 1][1]
			if(cg[i + 1][0] == 2):
				self.rpdelcnt += cg[i + 1][1]
				self.rpreflen += self.rpdelcnt

			if(debugDetailedMode == 1): 
				print("specialcase: rp refbaselen==self.rpreflen")
		else:
			'''bug:  for this special case.    fplen=20,    cigar:  15M8D 80M => reflen=13,seqlen=15
			'''
			if cg[i][0] == 2:
				self.rpreflen = refbaselen

			if(debugDetailedMode == 1): 
				print("specialcase: rp refbaselen > self.rpreflen and the count of deletion base is large")

		self.rpseqlen = self.rpreflen + self.rpinscnt - self.rpdelcnt

	def update(self, cigar):
		self.updateFPSeqLen(cigar)
		self.updateRPSeqLen(cigar)

		if(debugDetailedMode == 1): 
			self.displayprimerinfo()

	def displayprimerinfo(self):
		print('primer info fp(reflen,seqlen,inscnt,delcnt):(', self.fpreflen, ',', self.fpseqlen, ',', self.fpinscnt, ',', self.fpdelcnt, ')'),
		print('primer info rp(reflen,seqlen,inscnt,delcnt):(', self.rpreflen, ',', self.rpseqlen, ',', self.rpinscnt, ',', self.rpdelcnt, ')\n')

class AmpliconInfo:
	def __init__(self):
		self.no = ''
		self.fp = ''
		self.rp = ''
		self.gene = ''
		self.chrom = ''
		self.amps = 0
		self.inss = 0
		self.inse = 0
		self.ampe = 0
		self.fragmentlen = 0
		self.inslen = 0

	def update(self, no, fp, rp, gene, chrom, amps, inss, inse, ampe):
		self.no = no
		self.fp = fp
		self.rp = rp
		self.gene = gene
		self.chrom = chrom
		self.amps = int(amps)
		self.inss = int(inss)
		self.inse = int(inse)
		self.ampe = int(ampe)
		
		if(self.ampe > self.amps):
			self.fragmentlen = self.ampe - self.amps + 1
		else:
			self.fragmentlen = -1
			print("Warning:fragmentlen equals -1!")

		if(self.inse > self.inss):
			self.inslen = self.inse - self.inss + 1
		else:
			self.inslen = -1
			print("ERROR: Amplicon File is wrong!")

	def displayampliconinfo(self):
		print('amplicon info(no,fp,rp,gene,chrom,ampls,inss,inse,ample,fragment,lenfp,lenrp,leninsertion):(',self.no,',',self.fp,',',self.rp,',',self.gene,',',self.chrom,',',self.amps,',',self.inss,',',self.inse,',',self.ampe,',',self.fragmentlen,',',len(self.fp),',',len(self.rp),',',self.inslen,')')

def GetAmpliconsInfo(ampliconfile):
	ampll = list()
#	ampld = {}
	pf = open(ampliconfile)
	idx = 0

	for line in pf:
		line = line.strip('\n')
		frprimer = line.split(',')
		if(len(frprimer) == ampliconfilecol):
			idx += 1
			ampl = AmpliconInfo()
			ampl.update(frprimer[0], frprimer[1],frprimer[2], frprimer[3], frprimer[4], frprimer[5], frprimer[6], frprimer[7], frprimer[8])
			ampll.append(ampl)
#			ampld[idx] = ampl
	pf.close()

	return ampll
#	return ampld

## Process softclip situation, if the head of a read has softclip or the end has softclip.
def GetSeqandQualByDelSoftclip(cigar, seq, qual):
	presftclplen = backsftclplen = 0
	precgsftclplen = backcgsftclplen = 0
	seqnew = seq
	qualnew = qual
	cigarnew = cigar

	try:
		## ex: 10S130M
	        	## read.seq[10:150 - 10]
		if(cigar[0][0] == 4):
			presftclplen = cigar[0][1]
			precgsftclplen = 1
		## ex: 1. 130M20S, read.seq[0:150-20]   2. 10S120M20S,read.seq[10:150-20]
		## using 'precliplen - rdseqstart' to process ex.2. If don`t consider ex.2, rdseqstart = 0 is OK.
		if(cigar[len(cigar) - 1][0] == 4):
			backsftclplen = cigar[len(cigar) - 1][1]
			backcgsftclplen = 1
	except:
		print("Warning:GetSeqandQualByDelSoftclip except!")

	seqnew = seq[(0 + presftclplen):(len(seq) -backsftclplen)]
	qualnew = qual[(0 + presftclplen):(len(seq) -backsftclplen)]
	cigarnew = cigar[( 0 + precgsftclplen):(len(cigar) - backcgsftclplen)]

	if(debugPreProcessofSoftclip == 1):
		print("########pre softclip len:" + str(presftclplen) + ";back softclip len:" + str(backsftclplen) + ";GetSeqandQualByDelSoftclip")
		print("########old cigar:" + str(cigar) + ";old cigar:" + str(cigarnew))
		print("########odd seq:" + seq + ";new seq:" + seqnew)
		print("########odd qual:" + ''.join(chr(i+33) for i in qual) + ";new qual:" + ''.join(chr(i+33) for i in qualnew))

	return [cigarnew, seqnew, qualnew]

## compare primer with read sequence
def ComparePrimerwithSeq(reada, readb, fp, rp, fragmentlen, errCnt):
	rdaseq = reada.seq
	rdbseq = readb.seq
	isprimermatch = 'n'

	#### the shortest effecitive length of one read is 20.That means for one read, if we exclude the fp and rp,then the length must be more than 20.
	## for case: fp and rp are both in one read. ce chuan le
	if fragmentlen <= readseqlength and (len(rdaseq) < (len(fp) + len(rp) + 20) or len(rdbseq) < (len(fp) + len(rp) + 20)):
		return isprimermatch
	## for case: fp is in read+,part of rp or no rp in read+.
	if fragmentlen > readseqlength:
		if (reada.strand == '+' and len(rdaseq) < (len(fp) + 20)) or (readb.strand == '-' and (len(rdbseq) < len(rp) + 20)):
			return isprimermatch
		if (readb.strand == '+' and len(rdbseq) < (len(fp) + 20)) or (reada.strand == '-' and (len(rdaseq) < len(rp) + 20)):
			return isprimermatch

	## reada is + strand and readb is - strand. fp match reada and rp match readb.
	if(reada.strand == '+'):
		if regex.search(r'^(' + fp + '){e<=' + errCnt+ '}', rdaseq) and regex.search(r'^('+rp+'){e<='+errCnt+'}', DNA_revcplm(rdbseq)):
			isprimermatch = 'y'
	## reada is - strand and readb is + strand. fp match reada`s reverse seq and fp match readb.
	elif(reada.strand == '-'):
		if regex.search(r'^('+rp+'){e<='+errCnt+'}', DNA_revcplm(rdaseq)) and regex.search(r'^('+fp+'){e<='+errCnt+'}', rdbseq):
			isprimermatch = 'y'

	return isprimermatch

def GetReadCutPos(seqinfo, lenofread):
	## init
	rdcutstartpos = 0
	rdcutendpos = lenofread

	rdcutstartpos = seqinfo.fpseqlen
	rdcutendpos = lenofread - seqinfo.rpseqlen

	if(debugDetailedMode == 1):
		print('###### read cut pos:(start,end)=(', rdcutstartpos, ',', rdcutendpos, ')')

	return [rdcutstartpos, rdcutendpos]

def ModifyReadseqAndReadqual(oldrdseq, oldrdqual, cutstartpos, cutendpos):
	newrdseq = newrdqual = ''

	newrdseq = oldrdseq[cutstartpos:cutendpos]
	newrdqual = oldrdqual[cutstartpos:cutendpos]

	return [newrdseq, newrdqual]

def ModifyMapStartPos(oldmapspos, seqinfo):
	## init
	newmapspos = oldmapspos
	cutfpreflen = seqinfo.fpreflen

	newmapspos = oldmapspos + cutfpreflen

	return newmapspos

def ModifyMapEndPos(oldmapepos, seqinfo):
	## init
	newmapepos = oldmapepos
	cutrpreflen = seqinfo.rpreflen

	newmapepos = oldmapepos - cutrpreflen

	if(debugDetailedMode == 1):
		print('modify read map end pos:mapendpos=', newmapepos)
	return newmapepos

'''
Parameter: rd1info: (oldmapspos, seqinfo) rd2info:(oldmapspos, seqinfo)
'''
def ModifyReadPairMapStartPos(rdainfo, rdbinfo):
	rdaoldmapspos = rdainfo[0]
	rdaseqinfo = rdainfo[1]
	rdboldmapspos = rdbinfo[0]
	rdbseqinfo = rdbinfo[1]

	rdanewmapspos = ModifyMapStartPos(rdaoldmapspos, rdaseqinfo)
	rdbnewmapspos = ModifyMapStartPos(rdboldmapspos, rdbseqinfo)

	rdanewnextmapspos = rdbnewmapspos
	rdbnewnextmapspos = rdanewmapspos

	if(debugDetailedMode == 1):
		print('modify readpair map start pos:rdanew(mapspos,nextmapspos)=(',rdanewmapspos,',',rdanewnextmapspos,'),' , 'rdbnew(mapspos,nextmapspos)=(', rdbnewmapspos, ',', rdbnewnextmapspos, ')')

	return [(rdanewmapspos, rdanewnextmapspos), (rdbnewmapspos, rdbnewnextmapspos)]

'''
rd1info: (newrd1.strand, newrd1.mapspos, newrd1.mapepos) rd2info:(newrd2.strand, newrd2.mapspos, newrd2.mapepos)
'''
def ModifyReadPairTlen(rdainfo, rdbinfo):
	rdastrand = rdainfo[0]
	rdamapspos = rdainfo[1]
	rdamapepos = rdainfo[2]
	rdbstrand = rdbinfo[0]
	rdbmapspos = rdbinfo[1]
	rdbmapepos = rdbinfo[2]

	if(rdastrand == '+'):
		rdatlen = rdbmapepos - rdamapspos + 1
		rdbtlen = 0 - rdatlen
	else:
		rdbtlen = rdamapepos - rdbmapspos + 1
		rdatlen = 0 - rdbtlen

	if(debugDetailedMode == 1):
		print('modify tlen:rdanew=',rdatlen, ',rdbnew=', rdbtlen, ',rdastrand=', rdastrand, ',rdbstrand=', rdbstrand)

	return [rdatlen, rdbtlen]

def CutCigar(oldcigar, cutlen):
	baselen = 0
	index = -1
	newcigar = oldcigar
	cg = ()
	cgnum = 0

	'''10M1I3M2D5M1I 80M1I10M1I3M2D5M
	   0      1  2   3   4  5 6
	'''
	for x in oldcigar:
		index += 1
		if x[0] == 0 or x[0] ==1 or x[0] == 2:
			baselen += x[1]
			if(baselen >= cutlen):
				cgnum = baselen - cutlen
				break

	newcigar = oldcigar[(index+1):len(oldcigar)]
	if(cgnum != 0):
		'''for this case: 50M1I80M(fplen=20) => 30M1I80M   cgnum=30
		'''
		cg = (x[0], cgnum)
		newcigar.insert(0, cg)

	return newcigar

def ModifyCigar(oldcigar, seqinfo):
	tempcigar = newcigar = oldcigar

	fpcutcglen =  seqinfo.fpseqlen + seqinfo.fpdelcnt
	rpcutcglen =  seqinfo.rpseqlen + seqinfo.rpdelcnt

	tempcigar = CutCigar(oldcigar, fpcutcglen)
	newcigartemp = CutCigar(tempcigar[::-1],  rpcutcglen)
	
	newcigar = newcigartemp[::-1]

	if(debugDetailedMode == 1):
		print('oldcigar:', oldcigar,',tempcigar:',tempcigar, ',newcigar:', newcigar, ',fpcutlen:', fpcutcglen, ',rpcutlen:',rpcutcglen)

	return newcigar

'''
when cut fp:  cigar=['10','A','5','^','T','G','80','T','^','T','G','10','A','5']
when cut rp: cigar=['5', 'A', '10', 'G', 'T', '^', 'T', '80', 'G', 'T', '^', '5', 'A', '10']
If rp has ^, then ^ position is different from fp. This case is needed to process specially.
'''
def CutMD(oldmdarray, cutlen):
	index = 0
	totalcnt = 0
	new1ststartnum = 0
	newmdarray = oldmdarray

	while(index < len(oldmdarray)):
		if(oldmdarray[index].isdigit()):
			totalcnt += int(oldmdarray[index])
			if(totalcnt >= cutlen):
				new1ststartnum = totalcnt - cutlen
				break
			else:
				index += 1
		elif(oldmdarray[index].isalpha()):
			totalcnt += 1
			if(totalcnt >= cutlen):
				new1ststartnum = 0
				break
			else:
				index += 1
		elif(oldmdarray[index] == '^'):
			index += 1
			## while for this case, it has lots of deletion bases.   20M^ACGAAT80M
			while(oldmdarray[index].isalpha()):
				totalcnt += 1
				if(totalcnt >= cutlen):
					new1ststartnum = 0
					break
				else:
					index += 1
		else:
			print('ERROR:remove part MD!')
			return (-1,-1)			

	newmdarray = oldmdarray[(index + 1):len(oldmdarray)]

	'''only for cutting RP, because if RP is reversed, then deletion will be different compared with FP.
	  when cut rp: cigar=['5', 'A', '10', 'G', 'T', '^', 'T', '80', 'G', 'T', '^', '5', 'A', '10']      cutlen=18, index=4
	'''
	if (index + 1) < len(oldmdarray) and oldmdarray[index + 1] == '^' and (index + 2) <= len(oldmdarray):
		newmdarray = oldmdarray[(index + 2):len(oldmdarray)]

	if(new1ststartnum != 0):
		newmdarray.insert(0, str(new1ststartnum))

	if(debugDetailedMode == 1):
		print('cut MD(cutlen,index,totalcnt,new1ststartnum):', cutlen, index, totalcnt, new1ststartnum)

	return newmdarray

def ModifyMD(oldmd, seqinfo):
	oldmdarray = oldmd
	newmdarray = tempmdarray = oldmdarray

	fpcutmdlen =  seqinfo.fpreflen
	rpcutmdlen =  seqinfo.rpreflen

	tempmdarray = CutMD(oldmd, fpcutmdlen)
	newmdarraytemp = CutMD(tempmdarray[::-1],  rpcutmdlen)
	newmdarray = newmdarraytemp[::-1]

	if(debugDetailedMode == 1):
		print('oldmd:', oldmdarray,',tempmdarray:',tempmdarray, ',newmdarray:', newmdarray, ',fpcutlen:', fpcutmdlen, ',rpcutlen:',rpcutmdlen)

	return newmdarray

def ModifyASXS(oldAS, oldXS, oldmapspos, oldmapepos, newreflen):
	oldreflen = oldmapepos - oldmapspos + 1
	newAS = oldAS
	newXS = oldXS

	if oldreflen != 0:
		newAS = (oldAS * newreflen)/oldreflen
		newXS = (oldXS * newreflen)/oldreflen
	else:
		print('ERROR: oldreflen is 0!')

	if debugDetailedMode == 1:
		print('oldAS, oldreflen, newreflen', oldAS, oldreflen, newreflen)

	return [newAS, newXS]

def ModifyPartRead(rd1sinfo, rd2sinfo, ampl):
	rd1 = rd1sinfo
	rd2 = rd2sinfo
	newrd1 = ReadSimpleInfo()
	newrd2 = ReadSimpleInfo()

	if (debugSimpleMode == 1):
		print("####################begin to modify read pair######################")
		print('#########older info:seq,cigar,mapspos,mapepos,tlen,MD,NM,AS,XS')
		print('#########reada:',rd1.seq,rd1.cigar, rd1.mapspos, rd1.mapepos,rd1.tlen,rd1.MD,rd1.NM,rd1.AS,rd1.XS)
		print('#########readb:',rd2.seq,rd2.cigar, rd2.mapspos, rd2.mapepos,rd2.tlen,rd2.MD,rd2.NM,rd2.AS,rd2.XS)

	## init and update
	rd1seqinfo = SeqInfoInRead(rd1.mapspos, rd1.mapepos, ampl.inss, ampl.inse)
	rd2seqinfo = SeqInfoInRead(rd2.mapspos, rd2.mapepos, ampl.inss, ampl.inse)
	rd1seqinfo.update(rd1.cigar)
	rd2seqinfo.update(rd2.cigar)

	## modify sequence and quality
	rd1cutpos = GetReadCutPos(rd1seqinfo, len(rd1.seq))
	rd2cutpos = GetReadCutPos(rd2seqinfo, len(rd2.seq))

	if rd1cutpos[0] < 0 or rd1cutpos[1] < 0 or rd2cutpos[0] < 0 or rd2cutpos[1] < 0:
		return [newrd1, newrd2]

	newrd1.seq = rd1.seq[rd1cutpos[0]:rd1cutpos[1]]
	newrd1.qual = rd1.qual[rd1cutpos[0]:rd1cutpos[1]]
	newrd2.seq = rd2.seq[rd2cutpos[0]:rd2cutpos[1]]
	newrd2.qual = rd2.qual[rd2cutpos[0]:rd2cutpos[1]]

	## modify map-start-pos
	rd1rd2pos = ModifyReadPairMapStartPos((rd1.mapspos, rd1seqinfo), (rd2.mapspos, rd2seqinfo))
	newrd1.mapspos = rd1rd2pos[0][0]
	newrd1.nextmapspos = rd1rd2pos[0][1]

	newrd2.mapspos = rd1rd2pos[1][0]
	newrd2.nextmapspos = rd1rd2pos[1][1]

	## modify map-end-pos
	newrd1.mapepos = ModifyMapEndPos(rd1.mapepos, rd1seqinfo)
	newrd2.mapepos = ModifyMapEndPos(rd2.mapepos, rd2seqinfo)

	## modify tlen	
	newrd1.strand = rd1.strand
	newrd2.strand = rd2.strand
	rd1rd2tlen = ModifyReadPairTlen((newrd1.strand, newrd1.mapspos, newrd1.mapepos), (newrd2.strand, newrd2.mapspos, newrd2.mapepos))
	newrd1.tlen = rd1rd2tlen[0]
	newrd2.tlen = rd1rd2tlen[1]

	## modify cigar
	rd1cigar = ModifyCigar(rd1.cigar, rd1seqinfo)
	newrd1.cigar = rd1cigar

	rd2cigar = ModifyCigar(rd2.cigar, rd2seqinfo)
	newrd2.cigar = rd2cigar
	
	## modify MD
	newrd1.MD = ModifyMD(rd1.MD, rd1seqinfo)
	newrd2.MD = ModifyMD(rd2.MD, rd2seqinfo)

	## modify NM	
	newrd1.modifyNM()
	newrd2.modifyNM()

	## modify AS,XS 
	rd1ASXS = ModifyASXS(rd1.AS, rd1.XS, rd1.mapspos, rd1.mapepos, len(newrd1.seq))
	rd2ASXS = ModifyASXS(rd2.AS, rd2.XS, rd2.mapspos, rd2.mapepos, len(newrd2.seq))
	newrd1.AS = rd1ASXS[0]
	newrd1.XS = rd1ASXS[1]
	newrd2.AS = rd2ASXS[0]
	newrd2.XS = rd2ASXS[1]

	ErrorCheck1(rd1.qname, newrd1.seq, newrd1.cigar, newrd1.MD)
	ErrorCheck1(rd2.qname, newrd2.seq, newrd2.cigar, newrd2.MD)

	if (debugSimpleMode == 1):
		print('############new info:seq,cigar,mapspos,nextmapspos,tlen,MD,NM,AS,XS')
		print('#########reada:',newrd1.seq,newrd1.cigar, newrd1.mapspos, newrd1.nextmapspos,newrd1.tlen,newrd1.MD,newrd1.NM,newrd1.AS,newrd1.XS)
		print('#########readb:',newrd2.seq,newrd2.cigar, newrd2.mapspos, newrd2.nextmapspos,newrd2.tlen,newrd2.MD,newrd2.NM,newrd2.AS,newrd2.XS,'\n')

	return [newrd1, newrd2]

def ModifyRead(read, newreadsinfo):
	read.query_sequence = newreadsinfo.seq
	read.query_qualities = newreadsinfo.qual
	read.reference_start = newreadsinfo.mapspos - 1
	read.next_reference_start = newreadsinfo.nextmapspos - 1
	read.cigartuples = newreadsinfo.cigar
	read.template_length = newreadsinfo.tlen

	if True == read.has_tag('RG'):
		read.set_tags([('NM', newreadsinfo.NM), ('MD', ''.join(newreadsinfo.MD)), ('AS', str(newreadsinfo.AS)), ('XS', str(newreadsinfo.XS)), ('RG', read.get_tag('RG'))])
	else:
		read.set_tags([('NM', newreadsinfo.NM), ('MD', ''.join(newreadsinfo.MD)), ('AS', str(newreadsinfo.AS)), ('XS', str(newreadsinfo.XS))])

def ErrorCheck1(qname, seq, cigar, md):
	seqlen = len(seq)

	cg_match_ins_len = 0
	cg_match_del_len = 0
	for x in cigar:
		if x[0] == 0:
			cg_match_ins_len += x[1]
			cg_match_del_len += x[1]
		if x[0] == 1:
			cg_match_ins_len += x[1]
		if x[0] == 2:
			cg_match_del_len += x[1]

	mdlen = 0
	for x in md:
		if x.isdigit():
			mdlen += int(x)
		if x.isalpha():
			mdlen += 1
		if x == '^':
			mdlen += 0

	if cigar[0][0] == 1:
		print('ERROR:CigarFirstIsInsertion')
	if md[0].isalpha():
		if debugDetailedMode == 1:
			print('MDFirstIsMismatch!')
	if md[0] == '^':
		print('ERROR:MDFirstIsDeletion!')
	
	if seqlen != cg_match_ins_len or cg_match_del_len != mdlen:
		print('ERROR:length check NG! qname,seqlen,cg_m_i_len,cg_m_d_len,mdlen', qname, seqlen, cg_match_ins_len, cg_match_del_len, mdlen)



def GetModifyReadfor1Ampl(bamf, outputdir, ampls, errCnt, outputq, lock):
	print("GetModifyReadfor1Ampl:amplno=", ampls[0].no)
	rsltl = []
	outbfname = outputdir + '/' +  os.path.split(bamf)[1] + '.' + ampls[0].no + '.cutp.part.bam'

	samf = pysam.AlignmentFile(bamf, "rb")
	outbamf = pysam.AlignmentFile(outbfname, "wb", template = samf)

	for ampl in ampls:
		eachaplcnt = 0
		readdict = {}
		for read in samf.fetch(ampl.chrom, ampl.amps, ampl.ampe):
			l = list()
			if read.mapping_quality and not(read.is_secondary) and not(read.is_supplementary):
				l = GetSeqandQualByDelSoftclip(read.cigartuples, read.query_sequence, read.query_qualities)
				read.cigartuples = l[0]
				read.query_sequence = l[1]
				read.query_qualities = l[2]
	
				## check whether this readname is processed
	#			if read.qname in rmdname:
	#				continue
	
				if read.qname not in readdict:
					readdict[read.qname] = read
				else:
					rda = ReadSimpleInfo()
					rdb = ReadSimpleInfo()
				
					rda.update(read)
					rdb.update(readdict[read.qname])
	
					ispm = ComparePrimerwithSeq(rda,  rdb, ampl.fp, ampl.rp, ampl.fragmentlen, errCnt)
					if(ispm == 'y'):
						if(debugSimpleMode == 1):
							ampl.displayampliconinfo()
							print('\nread pair is match!,read name:', read.qname)
						
						readpairinfo = list()
						readpairinfo = ModifyPartRead(rda, rdb, ampl)
	
						if len(readpairinfo[0].seq) != 0 and len(readpairinfo[1].seq) != 0:
							ModifyRead(read, readpairinfo[0])
							ModifyRead(readdict[read.qname], readpairinfo[1])
	
							outbamf.write(read)
							outbamf.write(readdict[read.qname])
	
	#						lock.acquire()
	#						try:
	#							print('lock is get!')
	#						finally:
	#							lock.release()
	
							eachaplcnt += 2
	
							del readdict[read.qname]
	
						#rdname[read.qname] = 1
			#else:
			#	outnotmatchbamf.write(read)
		rsltl.append((ampl,eachaplcnt))

	outputq.put((rsltl, outbfname))

	samf.close()
	outbamf.close()



def GetAmpliconsLoopIndex(length, times):
	stidx = endidx = []
	cntof1time = round(length/times)

	if cntof1time <= 0:
		cntof1time = 1

	stidx = list(range(0, length, cntof1time))
	endidx = list(range(cntof1time, length, cntof1time))
	endidx.append(length)

	return list(zip(stidx, endidx))


def main():
	parser = argparse.ArgumentParser(description='remove primer seq of the read')
	parser.add_argument('-a', help='amplicon file containing primer info. Format(9 columns): No,fp,rp,gene,chrom,ampl_amplstart,ampl_insertionstart,ampl_insertionend,ampl_amplend')
	parser.add_argument('-b', help='bam file to be processed, must be ended with .bam')
	parser.add_argument('-o', help='output file')
	parser.add_argument('-l', help='readseq length in the sequencing platform')
	parser.add_argument('-t', help='threads')
	parser.add_argument('-e', help='error count for comparing each primer')
	args = parser.parse_args()

	primerfile = args.a
	bamf = args.b
	output = args.o
	threads = args.t
	errCnt = args.e

	global readseqlength
	readseqlength = int(args.l)

	outpath = os.path.split(output)[0]
	outname = os.path.split(output)[1]
	tempdir = outpath + '/temp'
	MakeDir(outpath)
	MakeDir(tempdir)

	if 1 != IsBamFile(bamf):
		print("ERROR:please input bam file!")
		return

	print("processed bam file is : " + bamf)

	## get primer list
	amplsinfo = GetAmpliconsInfo(primerfile)
	## divide amplicons into 10 pieces
	ampldivideidx = GetAmpliconsLoopIndex(len(amplsinfo), 10)

	if(debugSimpleMode == 1):
		for x in ampldivideidx:
			print(x[0], x[1], amplsinfo[x[0]].no)

	outamplcntf = open(tempdir + '/' + outname.split('.')[0] + '.amplcnt', "w")

	## Attention: must use Manager.Queue(), or else the programe would not exit if using mp.Queue()
	mg = Manager()
	output = mg.Queue()
	lock = Lock()
	amplcnts = []
	outpartbfnames = []

	processes = [mp.Process(target=GetModifyReadfor1Ampl, args=(bamf, tempdir, amplsinfo[idx[0]:idx[1]], errCnt, output, lock)) for idx in ampldivideidx]

	for p in processes:
		p.start()

	for p in processes:
		partrslt = output.get()
		'''get amplicons` count list and output bamfile names
		'''
		for x in partrslt[0]:
			amplcnts.append(x)
		outpartbfnames.append(partrslt[1])

	for p in processes:
		p.join()

	## merge all part bams
	mgbfname = outpath + '/' + outname + '.rmp.bam'
	subprocess.call("samtools merge -f -h " + bamf + " " +  mgbfname + " " + ' '.join(outpartbfnames), shell=True)

#	if debugSimpleMode == 1:
#		print("samtools merge -f -h " + bamf + " " +  mgbfname + " " + ' '.join(outpartbfnames))

	## sort remove-primer bam files
	subprocess.call("samtools sort -@ " + threads + ' ' + mgbfname + ' -f ' + mgbfname + '.sort.bam', shell = True)

	if os.path.exists(mgbfname) and os.path.exists(mgbfname + '.sort.bam'):
		os.remove(mgbfname)
		os.rename(mgbfname + '.sort.bam', mgbfname)

	## index bam files
	subprocess.call("samtools index " + mgbfname, shell = True)

	## delete all part bam files
	for x in outpartbfnames:
		if os.path.exists(x):
			os.remove(x)

	## write amplicons` count to one file
	for x in amplcnts:
		print('amplcnts:', x[0].no, x[1])
		outamplcntf.write(x[0].no + ',' + x[0].fp + ',' + x[0].rp + ',' + x[0].chrom + ',' + str(x[1]) + ',' + '\n')

	outamplcntf.close()


if __name__ == '__main__':
        main()


