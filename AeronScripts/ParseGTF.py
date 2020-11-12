#1       havana  exon    13221   14409   .       +       .       gene_id "ENSG00000223972"; gene_version "5"; transcript_id "ENST00000456328"; transcript_version "2"; exon_number "3"; gene_name "DDX11L1"; gene_source "havana"; gene_biotype "transcribed_unprocessed_pseudogene"; transcript_name "DDX11L1-202"; transcript_source "havana"; transcript_biotype "processed_transcript"; exon_id "ENSE00002312635"; exon_version "1"; tag "basic"; transcript_support_level "1";


import os
import sys
import re
from collections import defaultdict
from gtfparse import read_gtf
from tqdm import trange

es  = defaultdict(list) #retrieve exon start give an exon id
ee  = defaultdict(list) #retrieve exon end given an exon id
eg  = defaultdict(list) #retrieve exon ids given an gene id
ese = defaultdict(list)	#retrieve sequences given an exon id
ec  = defaultdict(list) #retrieve chromosome information for an exon id
est = defaultdict(list) #retrieve strand Information for an exon id
en  = defaultdict(list) #retrieve exon number for an exon id
et  = defaultdict(list)


class ParseGTF():
	def __init__(self, gtf):
		df = read_gtf(gtf)
		exons_df = df.loc[df["feature"] == "exon",:]
		exons_df = exons_df.reset_index(drop = True)
		exons_df["exon_id"] = exons_df["gene_id"] + exons_df["strand"]
		for i in trange(exons_df.shape[0]):
			exn = exons_df.loc[i, "exon_id"]
			gn = exons_df.loc[i, "gene_id"]
			tr = exons_df.loc[i, "transcript_id"]
			enu = exons_df.loc[i, "exon_number"]
			es[exn].append(exons_df.loc[i, "start"])
			ee[exn].append(exons_df.loc[i, "end"])
			ec[exn].append(exons_df.loc[i, "seqname"])
			est[exn].append(exons_df.loc[i, "strand"])
			et[exn].append(tr)
			eg[gn].append(exn)
			eg[tr].append(exn)
			en[exn].append(enu)

	def getExons(self, key):
		return eg[key]
	
	def getStart(self, key):
		print(es[key])
		return int(es[key][0])

	def getEnd(self, key):
		return int(ee[key][0])

	def getChromosome(self, key):
		return ec[key]

	def getStrand(self, key):
		return est[key]

	def getExonNumber(self, key):
		return int(en[key][0])

	def getTranscript(self, key):
		return et[key][0]

	def getAllGenes(self):
		i=0
		genes = []
		for key in eg.keys():
			if not (key.find("ENSG")):
				genes.append(key)
				i=i+1
		return genes

	def getTranscriptPosition(self, key):
		allexons = eg[key]
		startpoint=0
		endpoint=0
		if(len(allexons)>0):
			startpoint=int(ee[allexons[0]][0])
			endpoint=0
			for e in allexons:
				en=int(ee[e][0])
				if(en>endpoint):
					endpoint=en
		return startpoint,endpoint

	def getAllTranscripts(self):
		i=0
		transcripts = []
		for key in eg.keys():
			if not (key.find("ENST")):
				transcripts.append(key)
				i=i+1
		return transcripts


