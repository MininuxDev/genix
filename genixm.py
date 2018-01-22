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
from genixmf import *


from Bio import SeqIO
import os
import sys


#######Gestion des menus et choix

def options(*opt):
	
	"""Gestion des menus et des choix de l'utilisateur. Verification de la cohérence du choix, et gestion de la sortie de l'application. L'option de sortie 'x' est automatiquement ajoutée
	
	Argument: 
		
		*opt: dénomination des options sans leur numéro d'indice qui est fixé par l'ordre des options (a partir de 1)
		
	Return:
		
		si le choix est cohérent, retourne le numéro de l'option choisie - type str
	"""
	
	t = """Options:\n\n"""
				
	for i, j in enumerate(opt):
		
		t += "	{0}: {1} \n".format(i + 1, j)
	
	t += "	x: Retour/Quitter"
				
	choix = [str(i) for i in list(range(1,len(opt)+1))] + ["x"]
	
	r = None
	c = 0
	
	while r not in choix:
		
		if c < 8:
		
			print(t, '\n')
	
			r = input('>>> ')
		
			c += 1
			
		else:
			
			r = "x" 
			
	if r == "x":
		
		byebye()

	else:
		
		return r
	
	
	
	
def byebye():
	
	f = input("""Entrez "Q" pour quitter l'application ou "R" pour retourner au menu principal: """)
		
	if f == 'Q':
		
		sys.exit()
		
	else:
		
		menu()
		
		
		
def selid(max = None, mini = None):
	
	"""selection d'id dans la liste de la session de travail
	
	Argument:
		
		max: int / None , nombre maximal d'id séléctionables
		mini: int / None, nombre minimal d'id selectionables
	"""
	
	print('\n')
	
	print_session()
	
	lr = input("Entrer les ID des séquences à séléctionner ou 'x' pour quitter l'application ou revenir au menu principal.\n>>> ")
		
	lr = lr.split(' ')
		
		
		
	if 'x' in lr: 
		
		byebye()
		
		

		
	lr2 = []

	lso = list(SeqIO.parse(workfile, "fasta"))
		
	for k in lr:
		
		if k not in [str(i) for i in range(len(lso))]:
		
			print('\n')
			print('ID {0} non valide !\n'.format(k))
			
		else: lr2.append(k)
			
	
	
	if mini != None and len(lr2) < mini:
			
			print('Vous devez entrer au moins {0} ID et au plus {1} ID valides\n\n'.format(mini, max))
			lr2 = selid(max,mini)
			
			
	
	elif max != None and len(lr2) > max:
			
			print('Vous devez entrer au moins {0} ID et au plus {1} ID valides\n\n'.format(mini, max))
			lr2 = selid(max,mini)	
	
			
	
	
	return lr2
			
							
			
			
		
		
		
		
		
#################

def licence():
	
	print('\n\n')
	
	txt = """
	======================================================================
	                                 Genix
	
	
	Copyright © 2017 Yann BOUYERON
	
	
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

	======================================================================"""

	print(txt)
	
	print('\n\n')
	
	menu()


def menu():
	
	print("""Genix - menu\n""")
	
	r = options('Gerez votre session de travail','Convertir des séquences','Créer des séquences','Obtenir des informations sur des séquences','Aligner des séquences')

	if r == "1": 
		session()
	
	elif r == "2": 
		convert_m()
	
	elif r == "3": 
		creat_m()
	
	elif r == "4": 
		info_m()
	
	elif r == "5": 
		aln_m()
	
	



def session():
	
	"""gestion des sequences"""
	
	if seqwork_is_empty_m():
		
		print_session()
						
		r = options('Rechercher et ajouter des sequences','Supprimer des sequences','Effacer la session de travail','Retour au menu principal')
		
		
	else: 
		
		print('\n')
		print("Votre session de travail est vide. Vous devez ajouter des séquences pour pouvoir commencer à travailler.\n")
		
		
		r = '1'
			
		

	if r == '1':
				
		print("Entrez les mots clés des séquences à ajouter.")
		n = input('>>> ')
		
		print('\n')
		
		print("Entrer le path du repertoire de recherche [/home/]")
		p = input('>>> ')
		
			
		if os.path.isdir(p) or os.path.isfile(p):
			
			path = p
			
		else:
			
			path = '/home/'
			
		
		
		n = n.split(' ')
	
		search(path, *n, add = True)
				
		session()
			
	elif r == '2': 
		
		print("Entrez les id des séquences à supprimer.")
		n = input('>>> ')
		
		print('\n')
		
		dels(*n)
		
		session()
			
	elif r == '3':
		
		with open(workfile,'wb') as sq:
			pass
			
		print('Votre session de travail a été vidée.')
		print('\n')
		
		session()
		
	
	else: menu()
	
	
	
	






def info_m():
	
	"""information sur les sequences"""

	print("""Genix - info\n""")
	
	lso = list(SeqIO.parse(workfile, "fasta"))
	
	x = selid(mini=1)
	
	for i in x:
		
		try: 
			
			j = int(i)

			if j in range(len(lso)):
				
				i = lso[j]

				a = testalpha(i.seq)
				
				print('\n')
				
				print('name: ', i.name)
				print('description: ', i.description)
				
				print('\n')

				
				if a == IUPAC.unambiguous_dna or a == IUPAC.unambiguous_rna or a == generic_nucleotide:

					print('alphabet: ', a, '\n')
					
					long = len(i.seq)
					
					print('longueur: ', str(long))
					print('\n')

					txt = show(i.seq)
					
					print(txt, '\n')

				elif a == IUPAC.protein or a == IUPAC.extended_protein:
					
					long = int(len(i.seq))
					
					options = input('Affichage de la séquence peptitidique à 1 ou 3 lettres par acides aminés: [3] ')
					
					
						
						
					
					print('alphabet: ', a, '\n')
					
						
					print('longueur: ', str(long))
					
					print('\n')
					
					if options == '1':
					
						txt = showp(i.seq)
						
					else:
						
						txt = showp(seq3(i.seq))
					
					print(txt, '\n')
			
			else:
				print("""{0} n'est pas un id de séquence !""".format(i))
			
		except ValueError:
			print("""{0} n'est pas un id de séquence !""".format(i))
			
	menu()
	


def creat_m():
	
	"""creation de sequences"""
	
	print("""Genix - creat\n""")
	
	name = input('sequence name: ')
	seq = input('sequence: ')
	description = input('description: ')
	
	if langage(seq) == 3:
	
		seq = seq1(seq)
	
	
		
	elif testalpha(seq) == None:
		
		print("Le langage de votre séquence n'est pas reconnu !")
		
		menu()
	
	else:
		
		seq = seq.upper()
	
	
	print('\n')
	
	out = input('Voulez vous sauvegarder cette séquence dans un fasta ? (y/n): ')
	
	if out == 'y': 
		out = True
	
	else: 
		out = False
	
	f = creat(name, seq, description, out)
	
	print('\n')
	
	print (f)
	
	menu()
	
	
	
	
def convert_m():
	
	"""menu generique pour la gestion des convertions de sequences"""
	
	if seqwork_is_empty_m() == False: 
		session()
	
	print("""Genix - convert\n""")
	
	r = options('Transcrire','Retro-transcrire','Traduire','Tables des codes génétiques')
			
	if r == "1": 
		transcrib_m()
	
	elif r == "2": 
		rtrans_m()
	
	elif r == "3": 
		translate_m()
		
	#affichage des tables du code génétique	
	elif r == "4": 
		
		d = CodonTable.unambiguous_dna_by_id
		r = CodonTable.unambiguous_rna_by_id
	
		for i in d:
				
			print(d[i].id, ':', d[i].names[0])
			
		print('\n')
		
		table = input("Entrez le numéro de la table à afficher: [1 (généirque table)] ")
			
		if table not in ['1','2','3','4','5','6','9','10','11','12','13','14','15','16','21','22','23','24','25','26']:
				
			table = 1
				
		else: 
				
			table = int(table)
			
		print('\n')
		
		print('DNA table: \n')
		
		print(d[table])
		
		print('\n')
		
		print('RNA table: \n')
		
		print(r[table])
		
		convert_m()
		
		
		
	
	








def transcrib_m():
	
	"""gestion des transcriptions"""
	
	if seqwork_is_empty_m() == False: 
		session()
		
	print("""Genix - transcription\n""")	
	
	lso = list(SeqIO.parse(workfile, "fasta"))
	
	r = selid(mini=1)
			
	for i in r:
		
		i = int(i)
		
		bnt_dna_seq = lso[i].seq
		
		if testalpha(bnt_dna_seq) == IUPAC.unambiguous_dna:
		
			dna_name = lso[i].name
			
			rna_name = dna_name + '.arn'
			
			des = "transcript of {0} [{1}]".format(dna_name, lso[i].description)
			
			rna_seq = bnt_dna_seq.transcribe()
				
			#creation d'un objet SeqRecord
			f = creat(rna_name, str(rna_seq), des)
			
			#affichage
			
			print('\n')
			
			print('Transcription de {0}'.format(dna_name))
			
			
			print('\n')
			
			print(show(rna_seq))
		
			print('\n')
			
			out = input('Voulez vous sauvegarder cette séquence dans un fasta ? (y/n): ')
		
			if out == 'y': 
				lso = list(SeqIO.parse(workfile, "fasta"))
				id = len(lso) - 1
				mkfas(os.getcwd(), id)
				
		else: 
			print("Vous ne pouvez pas transcrire une séquence de type {0}".format(testalpha(bnt_dna_seq)))
				
	menu()
	
	
				
				
				
				
								
				
				









def rtrans_m():
	
	"""gestion des retrotranscriptions"""
	
	if seqwork_is_empty_m() == False: 
		session()
	
	print("""Genix - retro-transcription\n""")	
		
	lso = list(SeqIO.parse(workfile, "fasta"))
					
	r = selid(mini=1)

	for i in r:
		
		i = int(i)
			
		rna_seq = lso[i].seq

		if testalpha(rna_seq) == IUPAC.unambiguous_rna:
		
			rna_name = lso[i].name
			
			cdna_name = rna_name + '.adn'
			
			des = "retro transcript of {0} [{1}]".format(rna_name, lso[i].description)
			
			cdna_seq = rna_seq.back_transcribe()
				
			#creation d'un objet SeqRecord
			f = creat(cdna_name, str(cdna_seq), des)
			
			#affichage
			
			print('\n')
			
			print('Rétro-transcription de {0}'.format(rna_name))
			
			
			print('\n')
			
			print(show(cdna_seq))
		
			print('\n')
			
			out = input('Voulez vous sauvegarder cette séquence dans un fasta ? (y/n): ')
		
			if out == 'y': 
				lso = list(SeqIO.parse(workfile, "fasta"))				
				id = len(lso) - 1
				mkfas(os.getcwd(), id)
				
		else: 	
			print("Vous ne pouvez pas retro transcrire une séquence de type {0}".format(testalpha(rna_seq)))
	
	menu()
				
	
	
	
	
def translate_m():
	
	"""gestion des traductions"""
	
	if seqwork_is_empty_m() == False: 
		session()
	
	print("""Genix - translation\n""")	
		
	lso = list(SeqIO.parse(workfile, "fasta"))
					
	r = selid(mini=1)

	for i in r:
		
		i = int(i)
		
		seq_to_translate = lso[i].seq
	
		
		if testalpha(seq_to_translate) == IUPAC.unambiguous_dna or testalpha(seq_to_translate) == IUPAC.unambiguous_rna or testalpha(seq_to_translate) == generic_nucleotide:
			
			table = input("Entrez le numéro de la table à utiliser: [1 (standard table)] ")
			
			if table not in ['1','2','3','4','5','6','9','10','11','12','13','14','15','16','21','22','23','24','25','26']:
				
				table = 1
				
			else: 
				
				table = int(table)
				
			p_name = lso[i].name + '.pro'
			
			des = "translation of {0} , use table {1}".format(lso[i].name, CodonTable.unambiguous_dna_by_id[table].names[0])
			
			#translation
			prot_seq = seq_to_translate.translate(table = table, to_stop=True)
				
				
			#option d'affichage 
				
			f = creat(p_name, str(prot_seq), des)
				
			#affichage
			
			print('\n')
			
			print('Traduction de {0} - {1}'.format(lso[i].name, CodonTable.unambiguous_dna_by_id[table].names[0]))
				
			print('\n')
			print(showp(seq3(prot_seq)))
		
			print('\n')
			
			out = input('Voulez vous sauvegarder cette séquence dans un fasta ? (y/n): ')
		
			if out == 'y': 
				lso = list(SeqIO.parse(workfile, "fasta"))				
				id = len(lso) - 1
				mkfas(os.getcwd(), id)
				
		else: 	
			print("Il n'est pas possible de traduire une séquence de type {0}".format(testalpha(seq_to_translate)))
	
	menu()
	
	
	
	
	
def aln_m():
	
	"""gestion des alignements de sequences"""
	
	if seqwork_is_empty_m() == False: 
		session()
		
	lso = list(SeqIO.parse(workfile, "fasta"))
		
	print("""Genix - alignement / comparaison\n""")	
					
	r = options('Alignement global par Needleman & Wunsch','Alignement local par Smith & Waterman','Alignement multiple par clustalw')
		
	
	
	
	
	#alignement needle ou water avec emboss		
	if r == "1" or r == "2":
		
		id = selid(max = 2, mini = 2)
		
		if compatibilite_aln(*id) == False:
			
			print('\n')
			print("Il n'est pas possible d'aligner une sequence nucléotidique avec une séquence peptidique ! \n")
			aln_m()
		
		gop = input("""pénalité de délétion [10]: """)
		gex = input("""pénalité d'extension de délétion [0.5]: """)
		
		try:
			gop = float(gop)
			
		except ValueError: 
			gop = 10
		
		try:
			gex = float(gex)
			
		except ValueError:
			gex = 0.5
		
		
		#realisation de l'alignement, dans les deux cas (needle ou water), le fichier créé est le meme : emb.aln
		if r == "1":
			needle(*id, gop=gop, gex=gex)
			program = "Needleman & Wunsch"
		
		elif r == "2":	
			water(*id, gop=gop, gex=gex)
			program = "Smith & Waterman"
			
		
		#instenciation d'un objet aln
		aln = Emboss_aln_reader('emb.aln')
				
		#affichage de l'alignement
		print('\n')
		
		print("Choisissez votre format d'affichage\n\n")
		
		f = options('genix format','clustal format','emboss format','phylip format','anagene format')

		print('\n')

		#genix format
		if f == "1":
			
			sa = lso[int(id[0])].seq
			
			if testalpha(sa) == IUPAC.protein or testalpha(sa) == IUPAC.extended_protein: 
			
				x = showgenixp3(aln.aln, program=program, identity=aln.identity)
			
			else:
				
				x = showgenix(aln.aln, program=program, identity=aln.identity)
			
		#clustal format
		elif f == '2':
			x = aln.aln.format("clustal")
			
			
		#emboss format
		elif f == "3":
			x = aln.read_aln
			
			
			
		#phylip format
		elif f == "4":
			x = aln.aln.format("phylip")
			
		#anagene format
		elif f == "5":
			x = showanag(aln.aln, program=program, identity=aln.identity)
			
		print(x)
			
				
		menu()
		
		
		
		
	#alignement clustalw
	elif r == '3':
		
		id = selid(mini = 2)
		
		if compatibilite_aln(*id) == False:
			
			print('\n')
			print("Il n'est pas possible d'aligner des sequences nucléotidiques avec des séquences peptidiques ! \n")
			aln_m()
		
		#realisation de l'alignement
		align, stdout = clustal(*id)
		
		score_matrice = matrix_pd(stdout)
		score_str = str(score_matrice)
		score = score_str.split('\n')
		score = "\n        ".join(score)
		
		#affichage de l'alignement
		print('\n')
		
		print("Choisissez votre format d'affichage\n\n")
		
		f = options('Genix format','Clustal format','Phylip format','Anagene format')

		print('\n')
		
		
		#genix format
		if f == "1":
			
			sa = lso[int(id[0])].seq
			
			if testalpha(sa) == IUPAC.protein or testalpha(sa) == IUPAC.extended_protein: 
			
				x = showgenixp3(align, program="clustalw", identity=score)
			
			else:
			
				x = showgenix(align, program="clustalw", identity=score)
			
			
		
		#clustal format
		elif f == "2":
			x = align.format('clustal')
			
			
		#phylip format
		elif f == "3":
			x = align.format('phylip')
		
		#anagene format			
		elif f == "4":
			x = showanag(align, program="clustalw", identity=score)
			
		
		print(x)
		
		menu()
		
	
		
