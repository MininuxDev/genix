#!/usr/bin/python3
# -*- coding: utf-8 -*-
#
# Genix - Outil pédagogique d'analyse génétique et phylogénétique
#
#
# Copyright (c) 2017 Yann BOUYERON
#
#
# licensed under GNU GPL version 3 (or later)
#
#
#
# This file is part of Genix.
#
# Genix is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Genix is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Genix.  If not, see <http://www.gnu.org/licenses/>.



from genixf import *


from Bio.Emboss.Applications import NeedleCommandline
from Bio.Emboss.Applications import WaterCommandline
from Bio import SeqIO
from Bio.Align.Applications import ClustalwCommandline
from Bio import AlignIO
import os


#Ensemble de fonctions dédiées à genix menu

workfile = os.getcwd() + "/.seqwork"

		

def print_session():
	
	"""affiche le contenu du workfile"""
	
	lso = list(SeqIO.parse(workfile, "fasta"))
		
	print('\n')
	print("""Genix  -  session""")
	
	print('\n')
	print(("Contenu de votre session de travail \n"))

	for i,j in enumerate(lso):

		print(i, ' - Nom: ', j.name)
		print('Description: ', j.description, '\n')
		
		
		



def seqwork_is_empty_m2():
	
	"""determine si le workfile est plein
	return True si il est plein, et False si il est vide"""
	
	if os.path.isfile(workfile):
	
		lso = list(SeqIO.parse(workfile, "fasta"))
		
		if len(lso) > 0:
			return True
			
		else: 
			return False
	
	else: 
		return False


def seqwork_is_empty_m():
	
	"""determine si le workfile est plein
	return True si il est plein, et False si il est vide"""
	
	try: 
		
		lso = list(SeqIO.parse(workfile, "fasta"))
		
	except FileNotFoundError:
		
		return False
			
	else:
	
		if len(lso) > 0:
			return True
			
		else: 
			return False
	
	



def search(path, *keys, add = False):

	"""
	recherche et ajoute des seqrecord dans le workfile
	
	Arguments:
		path: str : path du repertoire de recherche
		keys: str : mots clés de recherche
		add: boléen: si True, chaque seqrecord recherché est ajouté dans workfile ; si False, workfile est ecrasé
	Return:
		une liste des seqrec ajoutés au fasta temporaire
	"""
	
	lso = []
	
	#si path ne se termine pas par un / on le rajoute de maniière à obtenir des path complets depuis la racine
	if path[len(path)-1] != '/':
		path = path + '/'
	
	#si la tuple de mot_cle est vide, on créé le mot_cle '' qui correspond à n'importe quel caractère
	if keys == ():
		keys = ('',)
	
	for i in recurliste(path):
	
		name = i.split('/')
		name = name[len(name)-1]
		
		if os.path.isfile(i) and '.fas' in name:
			for j in keys:
				if j in name:
					so = list(SeqIO.parse(i, "fasta"))
					lso = lso + so
					break  


	#elimination des doublons
	list_SeqRec_open = []
	
	list_name =[]
	
	for i in lso:
		
		if i.name not in list_name:
			list_SeqRec_open.append(i)
			list_name.append(i.name)

	

	#ecriture des seqrec dans workfile
	if add == False:
		
		#ecriture de la liste seqrec dans le workfile. si le fichier existe déjà il est ecrasé
		SeqIO.write(list_SeqRec_open, workfile, "fasta")
	
	if add == True:
		
		#ecriture de la liste seqrec à la suite du workfile
		with open(workfile, 'ab') as fil:
			
			for k in list_SeqRec_open:
				
				x = k.format('fasta')
				fil.write(x.encode())
	
	return list(SeqIO.parse(workfile, "fasta"))


def dels(*id):
	
	"""
	supprime des seqrecord du workfile d'apres leur id
	
	Argument: 
		id: index des seqrec à supprimer du workfile
	
	"""
	lso = list(SeqIO.parse(workfile, "fasta"))
	
	nlso = []
	
	for i, j in enumerate(lso):
		
		if str(i) not in id:
			
			nlso.append(j)
			
			
	SeqIO.write(nlso, workfile, "fasta")
	
	
	


def creat(seq_name, seq, description, out = False):
	
	"""créer un fasta mono séquence et ajoute le seqrec au fasta temporaire
	
	Arguments:
		seq_name: (str) nom de la séquence créée
		seq: (str) la séquence
		description: (str) déscription du seqrec
		
	retrun: 
		le seqrec créé
	"""
	
	type = testalpha(seq)
	
	#creation d'un objet SeqRecord
	f = SeqRecord(Seq(seq, type), id= seq_name, name = seq_name, description = description)
	
	x = f.format('fasta')
	
	with open(workfile, 'ab') as fil:
		
		fil.write(x.encode())
	
	if out == True:

		SeqIO.write(f, os.getcwd() + '/' + f.name + '.fas', "fasta")
		
	return f



	


def mkfas(path, *seq):
	
	"""
	crée des fasta mono séquence
	
	arguments: 
		*seq: str: id des seqrec selectionnés dans le workfile
		path: str: repertoire dans lequel seront enregistrés les fasta créés
	return:
		cette fonction ne rtourne rien, elle crée des fichiers fasta mono séquence pour chaque id de séquence entré en argument
		les fichiers sont enregistrés dans le path indiqué en argument
	"""
	
	lso = list(SeqIO.parse(workfile, "fasta"))
	
	if os.path.isdir(path):
		
		if path[len(path)-1] != '/':
			
			path = path + '/'
		
		for i in seq:
			
			if os.path.isfile(i) and '.fas' in i:
				
				lsf = list(SeqIO.parse(i, "fasta"))
				
				for j in lsf:
					
					SeqIO.write(j, path + j.name, "fasta")
				
			else:
			
				SeqIO.write(lso[int(i)], path + lso[int(i)].name + '.fas', "fasta")



def mkfasx(out, *seq):
	
	"""
	crée des fasta multi séquence
	
	arguments: 
		*seq: id des seqrec selectionnés dans le workfile
		out: str: path absolu du fasta mutli séquence qui sera créé
	return:
		cette fonction ne rtourne rien, elle crée un fichier fasta multi séquences contenant chaque séquences dont l'id a été  entré en argument
	"""
	
	lso = list(SeqIO.parse(workfile, "fasta"))
	
	rec = [lso[int(i)] for i in seq]
			
	SeqIO.write(rec, out, "fasta")
			



def info(*id):

	"""
	arguments:
		id: seqrec selectionnés dans le workfile
		return les info des seqrec selectionnés par leur id
	"""

	lso = list(SeqIO.parse(workfile, "fasta"))
	
	lsi = []
	
	for i in id:
	
		try:
			
			i = int(i)
			i = lso[i]
			
			if testalpha(i.seq) == IUPAC.unambiguous_dna or testalpha(i.seq) == IUPAC.unambiguous_rna or testalpha(i.seq) == generic_nucleotide:

				lsi.append((i.name,i.description,str(len(i.seq)),testalphabet(i.seq),i.seq))

			elif testalpha(i.seq) == IUPAC.protein or testalpha(i.seq) == IUPAC.extended_protein:

				lsi.append((i.name, i.description, str(int(len(i.seq)/3)), testalpha(i.seq), i.seq))


		except:
			
			print('\n')
			print("""{0} n'est pas un identifiant de séquence valable.""".format(i))
			print('\n')
			
			pass
			
	return lsi









def clustal(*id, out = 'comp.aln'):
	
	mkfasx(out,*id)

	cline = ClustalwCommandline("clustalw", infile=out, score='PERCENT')
	
	myStdout, myStderr = cline()
	
	align = AlignIO.read(out, "clustal")
	
	return align, myStdout


def needle(*id,gop=10,gex=0.5,out ='emb.aln'):

	"""Alignement global par la methode de Needleman"""
	
	lso = list(SeqIO.parse(workfile, "fasta"))
	
	mkfasx('seqa.fas', id[0])
	
	mkfasx('seqb.fas', *id[1:])
	
	needle_cline = NeedleCommandline(asequence='seqa.fas', bsequence='seqb.fas', gapopen=gop, gapextend=gex, outfile=out)

	stdout, stderr = needle_cline()

	os.remove('seqa.fas')
	os.remove('seqb.fas')

	if len(id) < 3:
		align = AlignIO.read(out, "emboss")
		return align



def water(*id,gop=10,gex=0.5,out ='emb.aln'):

	"""Alignement global par la methode de Needleman"""
	
	lso = list(SeqIO.parse(workfile, "fasta"))
	
	mkfasx('seqa.fas', id[0])
	
	mkfasx('seqb.fas', *id[1:])
	
	water_cline = WaterCommandline(asequence='seqa.fas', bsequence='seqb.fas', gapopen=gop, gapextend=gex, outfile=out)

	stdout, stderr = water_cline()

	os.remove('seqa.fas')
	os.remove('seqb.fas')

	if len(id) < 3:
		align = AlignIO.read(out, "emboss")
		return align



def compatibilite_aln(*id):
	
	"""test la compatibilité des séquences à aligner"""
	
	lso = list(SeqIO.parse(workfile, "fasta"))
	
	langlist = [testalpha(lso[int(i)].seq) for i in id]
	
	if (IUPAC.protein in langlist or IUPAC.extended_protein in langlist) and (IUPAC.unambiguous_dna in langlist or IUPAC.unambiguous_rna in langlist or generic_nucleotide in langlist):
		
		return False
		
	else:
		
		return True


 
