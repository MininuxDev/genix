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


import numpy as np
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC, generic_nucleotide, generic_protein, SingleLetterAlphabet
from Bio.Data import CodonTable
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import GC, seq1, seq3
from Bio import SeqIO
from io import StringIO
from Bio import AlignIO
from Bio import Phylo
import os
import fnmatch
from xhtml2pdf import pisa
from xhtml2pdf.default import DEFAULT_CSS
from pandas import DataFrame
import unicodedata




def recurliste(path):
	
	"""
	path est une chaine de caractere
	retourne une liste recursive des dossiers et fichiers contenus dans path.
	
	c'est une alternative a glob.glob('**',recursive=True) qui ne fonctionne pas en python < 3.5"""

	r = []
	for root, dir, files in os.walk(path):
		if root != path:
			r.append(root.replace(path+'/', ''))
		for f in fnmatch.filter(files, "*"):
			x = root+'/'+f
			r.append(x)
	return r


############ TEST DE L'ALPHABET 	###########
def canbedna(seq):
	
	"""can be a dna ?
	argument: 
		seq (str): la séquence à tester
	return True / False"""

	div = list(set(seq))
			
	for i in div:
		
		if i not in list(IUPAC.unambiguous_dna.letters) and i.upper() not in list(IUPAC.unambiguous_dna.letters):
			
			return False
			
	return True


def canberna(seq):
	
	"""can be an rna ?
	argument: 
		seq (str): la séquence à tester
	return True / False"""
	
	div = list(set(seq))
	
	for i in div:
		
		if i not in list(IUPAC.unambiguous_rna.letters) and i.upper() not in list(IUPAC.unambiguous_rna.letters):
			
			return False
			
	return True


def canbeexprot(seq):
	
	"""can be an extended protein ?
	argument: 
		seq (str): la séquence à tester
	return True / False"""
	
	div = list(set(seq))
	
	for i in div:
		
		if i not in list(IUPAC.extended_protein.letters) + ['*']:
			
			return False
			
	return True


def canbeprot(seq):
	
	"""can be a protein ?
	argument: 
		seq (str): la séquence à tester
	return True / False"""
	
	div = list(set(seq))
	
	for i in div:
		
		if i not in list(IUPAC.protein.letters) + ['*']:
			
			return False
			
	return True






def testalpha(seq):
	
	"""Test l'alphabet d'une sequence
	
	Argument:
		seq: (str) , la sequence à tester
	
	Retourne:
		generic_nucleotide si il y'a un doute entre adn ou arn
		IUPAC.unambigous_dna
		IUPAC.unambigous_rna
		IUPAC.protein
		IUPAC.extended_protein
		ou None si la sequence n'est pas reconnue
		
	"""
	
	if canbedna(seq) and canberna(seq):
	
		return generic_nucleotide
	
	elif canbedna(seq):
		
		return IUPAC.unambiguous_dna
		
	elif canberna(seq):
		
		return IUPAC.unambiguous_rna
		
	
	elif canbeprot(seq):
		
		return IUPAC.protein
		
	elif canbeexprot(seq):
		
		return IUPAC.extended_protein
		
	else:
		
		return None




###############################




def langage(seq):
	
	"""determine le mode d'ecriture a 1 ou 3 lettres des séquences peptidiques
	
	Argument:
		
		seq: la sequence à tester (str) ou (objet Seq)
	
	Retourne:
		1 (int) si langage à une lettre
		3 (int) si langage à trois lettres
		ou None si la séquence n'est pas reconnue comme une sequence proteique
	"""
	
	seq = str(seq)

	alpha = testalpha(seq)

	if seq.isupper() and (alpha == IUPAC.protein or alpha == IUPAC.extended_protein):

		return 1

	elif seq.isupper() == False and seq.islower() == False and (testalpha(seq1(seq)) == IUPAC.protein or testalpha(seq1(seq)) == IUPAC.extended_protein):

		return 3

	else:

		return None






def array_matrix(myStdout):
	
	"""transforme le myStdout retourné par Clustalcommandeline avec l'option score, en une matrice au format ndarray"""

	y = myStdout.split('\n')
	
	legende =[]
	s = []
	
	for i in y:
		
		if 'Sequence ' in i and 'Multiple' not in i and 'Pearson' not in i: legende.append(i)
			
		elif 'Sequences ' in i: s.append(i)
			
			
			
	matrix = []
	
	for i in s:	
		z = i.split(' ')
		z = z[1] + z[len(z)-1]
		z = z[1:]
		z = z.replace(':', ' ')
		z = z.replace(')', ' ')
		z = z.split(' ')
		
		for i, j in enumerate(z): z[i] = int(j)
		
		matrix.append(z)
		
	
	
	m = np.zeros((len(legende)+1,len(legende)+1))
	
	for i, j in enumerate(m[0]): 
		m[0][i] = i
		
	for i, j in enumerate(m[None][0]): 
		m[i][0] = i
	
	for i in matrix:
		
		m[i[0]][i[1]]=i[2]
		m[i[1]][i[0]]=i[2]
		
	for i, j in enumerate(m):
		
		for k, l in enumerate(j):
			
			if l == 0:
				
				m[i][k] = 100
		
	m[0][0] = 0
	
	
		
	return legende, m	





def matrix(myStdout):
	
	"transforme le myStdout retourné par Clustalcommandeline avec l'option score, en une matrice sous forme d'une liste de liste (liste de ligne)"

	y = myStdout.split('\n')
	
	
	#extraction des legendes et des scores
	
	legende =[]
	
	s = []
	
	for i in y:
		
		if 'Sequence ' in i and 'Multiple' not in i and 'Pearson' not in i:
			legende.append(i)
			
		elif 'Sequences ' in i:
			s.append(i)
			
	m = []
	
	for i in s:
		
		z = i.split(' ')
		
		z = z[1] + z[len(z)-1]
		
		z = z[1:]
		
		z = z.replace(':', ' ')
		
		z = z.replace(')', ' ')
		
		z = z.split(' ')
		
		for i, j in enumerate(z): 
			
			z[i] = int(j)
		
		m.append(z)
		
	
	#construction d'une matrice vide
	matrix = [[]]
	
	for l in legende:
	
		matrix.append([])
	
	for i in matrix:
		
		i.append(' ')
		
		for l in legende:
			
			i.append(' ')
	
	
	
	#remplissage des titres
	
	for i, j in enumerate(legende):
		
		#redecoupage du nom
		j = j.replace('Sequence ', '')
		j = j[:len(j)-12]
	
	
		matrix[0][i+1] = j
		
		matrix[i+1][0] = j
		
		
		
	#remplissage des scores
	
	for i in m:
		
		matrix[i[0]][i[1]] = i[2]
		matrix[i[1]][i[0]] = i[2]
		
	for k in matrix:
		
		for i, j in enumerate(k):
			
			if j == ' ': k[i] = 100
		
	matrix[0][0] = """Genix"""
	
	
	return matrix 
	
	
	
	
	
	
	
def matrix_pd(myStdout):
	
	"transforme le myStdout retourné par Clustalcommandeline avec l'option score, en une matrice pandas"

	y = myStdout.split('\n')
	
	
	# extraction des legendes et des scores
	
	l =[]
	
	s = []
	
	for i in y:
		
		if 'Sequence ' in i and 'Multiple' not in i and 'Pearson' not in i:
			
			t = i.split(' ')
			
			l.append(t[2])
			
			
			
			
		elif 'Sequences ' in i:
			
			s.append(i)
			
			
	# decoupage des str contenu dans score, de maniere a obtenir une liste de liste ou chaque sous liste contient 3 elements: id,id,score: exemple: [[1, 2, 93], [1, 3, 96], [2, 3, 90]] ... -> score de 93 entre seq1 et seq2
	
	m = []
	
	for i in s:
		
		z = i.split(' ')
		
		z = z[1] + z[len(z)-1]
		
		z = z[1:]
		
		z = z.replace(':', ' ')
		
		z = z.replace(')', ' ')
		
		z = z.split(' ')
		
		for i, j in enumerate(z): 
			
			z[i] = int(j)
		
		m.append(z)
		
		
	# modification des id des sequences dans m de maniere à compter a partir de 0... -> [[0, 1, 93], [0, 2, 96], [1, 2, 90]]
	
	for i in m:
		
		i[0] = i[0] - 1
		i[1] = i[1] - 1
	
	
	# creation d'une matrice pleine de 100... -> [[100, 100, 100], [100, 100, 100], [100, 100, 100]]
	
	mat = []
	
	for i in l:

		mat.append([100]*len(l))
	
	
	# remplissage de la matrice... -> [[100, 93, 96], [100, 100, 90], [100, 100, 100]]
	
	for k in m:

		mat[k[0]][k[1]] = k[2]
		mat[k[1]][k[0]] = k[2]
		


	# mise en forme avec pandas
	
	matrixpd = DataFrame(mat, index = l, columns = l)
	
	
	return matrixpd
	
	
	
	
def matrix2html(title, matrix, legende, auteur= '', text=''):
	
	"""
	transforme une matice au format liste de liste ou un ndarray en html
	"""

	td = "".join(["""<tr>{x}</tr>""".format(x="".join(["""<td>{j}</td>""".format(j=j) for j in i])) for i in matrix])
	
	
	table = """<table><caption>{0}</caption>{1}</table>""".format(legende, td)
	
	css = """
	table {
	border-collapse: collapse; 
	min-width: 80%;
	margin-left: 5%;
	margin-right: 5‰;
	font-family: times;
	}

	h3 {
	font-family: times;
	}


	td {
	border: 1px solid black;
	text-align: center;
	width:auto;
	}


	tr {
	height: 35px;	
	}

	corps {
	width: 100%;
	margin-left: 0%;
	margin-right: 0‰;
	font-family: times;
	}
	
	.auteur {
	text-align: right;
	margin-right: 15%;
	font-style: italic;
	font-size: 80%;	
	}
	"""
	
	
	html = """
	<!doctype html>
	<html>
		<head>
			<meta charset="utf-8">
			<title>Title</title>
			<style type="text/css">
			{0}
			</style>
		</head>
		<body>
			<h3>{1}</h3>
			<div class='corps'>
			{2} 
			<p class = 'auteur'>Realised with Genix <br> {3}</p>
			</div>
			<p>{4}</p>
		</body>
	</html>
	
	""".format(css, title, table, auteur, text)
	
	return html


	
def htmlpdf(html, pdf, css = ''):
	
	"""Transforme un html déjà structuré (doctype head html body) en pdf
	
	
	Arguments:
		
		html: path ou str du html déjà structuré
		pdf: path du fichier pdf qui sera créé
		css: (type str) c'est le css qui sera ajouté au css par default de pisa
		Attention: le css de la balise style du html est prioritaire sur le css entré en argument !
		
	"""
	
	#recuperation du html 
	if os.path.isfile(html):
		
		with open(html,'rb') as h:
		
			html = h.read().decode()
		
				
	if 'doctype' in html and 'head' in html and 'html' in html:	
				
		#ecriture du pdf
		with open(pdf,'wb') as p:
			
			pisa.CreatePDF(html.encode(), p, default_css = DEFAULT_CSS + css)

	









	
	
	
	





	


	
	





def show(seq, start = 0, stop = None, width = 60):
	
	"""présentation d'une séquence nucléotique decouppée en tronçons avec regle graduée.
	start et stop doivent etre des entiers (int); width doit etre un multiple de 10.
	retourne une str
	"""
	
	seq = str(seq)

	if width/10 != int(width/10):
		return "Value Error: le parametre width doit etre un multiple de 10"

	if stop == None: 
		stop = len(seq)

	elif type(stop) == type(int()): 
		pass

	else: 
		return "Value Error: le parametre stop doit etre un entier"


	
	if type(start) != type(int()):
		return "Value Error: le parametre start doit etre un entier"

	d = start
	
	f = d + width
	
	

	txt = ''

	tir = "----:----|" * int(width/10)

	while d <= stop:
		
		if f > stop:
			
			tir = tir[:len(seq[d:stop])]
		
		
		
		num = "          "
		
		for k in np.arange(10, width + 10, 10):
				
			num = num + """{0}{1}""".format("         " if k == 10 else "        " if  d+k <= 100 else "       " if d+k <= 1000 else "      ", d + k if d+k <= stop else "")
		
		
		if f <= stop:
			
			txt = txt + "          " + seq[d:f] + "\n" + "          " + tir + "\n" + num + "\n\n"
		
		else:
			
			txt = txt + "          " + seq[d:stop] + "\n" + "          " + tir + "\n" + num + "\n\n"
		
		d, f = d + width, f + width

	return txt
	
	
	
	



		








def showp(seq, start = 0, stop = None, width = None):
	
	"""présentation d'une séquence peptidique decouppée en tronçons avec regle graduée.
	start et stop doivent etre des entiers (int); width doit etre un multiple de 10.
	retourne une str"""
	
	seq = str(seq)
	
	if langage(seq) == 1:
		
		if width == None: width = 60
	

		if width/10 != int(width/10):
			return "Value Error: le parametre width doit etre un multiple de 10"
	
		if stop == None: 
			stop = len(seq)
	
		elif type(stop) == type(int()): 
			pass
	
		else: 
			return "Value Error: le parametre stop doit etre un entier"
	
	
		
		if type(start) != type(int()):
			return "Value Error: le parametre start doit etre un entier"
	
		
			
		d = start
		
		f = d + width
	
		txt = ''	
		
		
		
	
	
		tir = "----:----|" * int(width/10)
	
		while d <= stop:
			
			if f > stop:
				
				tir = tir[:len(seq[d:stop])]
			
			
			
			num = "          "
			
			for k in np.arange(10, width + 10, 10):
					
				num = num + """{0}{1}""".format("         " if k == 10 else "        " if  d+k <= 100 else "       " , d + k if d+k <= stop else "")
			
			
			if f <= stop:
				
				txt = txt + "          " + seq[d:f] + "\n" + "          " + tir + "\n" + num + "\n\n"
			
			else:
				
				txt = txt + "          " + seq[d:stop] + "\n" + "          " + tir + "\n" + num + "\n\n"
			
			d, f = d + width, f + width
	
	
	
		return txt
		
		
		
	elif langage(seq) == 3:
		
		if width == None: width = 30
		
		if width/10 != int(width/10):
			return "Value Error: le parametre width doit etre un multiple de 10"
	
		if stop == None: 
			stop = int(len(seq)/3)
	
		elif type(stop) == type(int()): 
			pass
	
		else: 
			return "Value Error: le parametre stop doit etre un entier"
	
	
		
		if type(start) != type(int()):
			return "Value Error: le parametre start doit etre un entier"
	
		
			
		d = start
		
		f = d + width
	
		txt = ''	
		
		
		
		
		#decoupage de la séquence
		w = 0
		h = 3
		
		z = []
		
		while h <= len(seq):
			
			z.append(seq[w:h])
			
			w, h = h, h + 3
		
		
		
		
		
		
		tir = " -  -  -  -  :  -  -  -  -  | " * int(width/10)
	
		while d <= stop:
			
			if f > stop:
				
				tir = tir[:len(''.join(z[d:stop]))]
			
			
			
			num = "          "
			
			for k in np.arange(10, width + 10, 10):
					
				num = num + """{0}{1}""".format("                            " if  d+k <= 100 or k == 10 else "                           ", d + k if d+k <= stop else "")
			
			
			if f <= stop:
				
				txt = txt + "          " + ''.join(z[d:f]) + "\n" + "          " + tir + "\n" + num + "\n\n"
			
			else:
				
				txt = txt + "          " + ''.join(z[d:stop]) + "\n" + "          " + tir + "\n" + num + "\n\n"
			
			d, f = d + width, f + width
	
	
	
		return txt
		












def showanag(align, start = 0, stop = None, width = 60, program=None, identity=None):
	
	
	"""présentation d'un alignement nucléotique decouppée en tronçons avec regle graduée , format anagene-like.
	start et stop doivent etre des entiers; width doit etre un multiple de 10
	
	Arguments:
		align: l'alignement (str) retourné par la fonction AlignIO.read de biopython: exemple : align = AlignIO.read(path, "clustal")
		
		start: (int) index du premier nucléotide représenté
		stop: (int) index du dernier nucléotide représenté
		width: (int) largeur des lignes; width doit etre un multiple de 10
		program: (str) nom du programme ayant réaliser l'alignement (facultatif)
		identity: (str, int ou float) % d'identité ou matrice de score 
	"""
	
	
	s = list(align)
	
	#sequence de reférence 
	sa = str(s[0].seq)
		
	#autres sequence de l'alignement
	ss = [str(i.seq) for i in s[1:]]



	#verification des options

	if width/10 != int(width/10):
		return "Value Error: le parametre width doit etre un multiple de 10"

	if stop == None: 
		stop = len(sa)

	elif type(stop) == type(int()) and stop <= len(sa): 
		pass

	else: 
		return "Value Error: le parametre stop doit etre un entier inférieur ou égal à la longueur de l'alignement"

	if type(start) != type(int()):
		return "Value Error: le parametre start doit etre un entier"

	
	#modification de sa: remplacement des - par des _
	san = ""
	
	for i in sa:
		if i == '-': san += "_"
		else: san += i

	#modif des sequences secondaires de ss, création des symboles *
	
	sss = []
	
	for sb in ss:
		
		sbn = ""
	
		for a , b in zip(sa,sb):
			
			if b == "-": 
				
				sbn += "_"
				
				
			elif b == a: 
				
				sbn += "-"
				
				
			elif b != a: 
				
				sbn += b
				
		sss.append(sbn)





	#récupération des noms des séquences
	
	#nom de la séquence de référence
	nsa = str(s[0].id)
	
	#liste des noms des autres séquences de l'alignement
	nnn = [str(i.id) for i in s[1:]]
	
	
	#recherche du nom de séquence le plus long et gestion des espaces entre nom et séquences
	lt = nnn
	lt.append(nsa)
	
	ltlen = []
	for i in lt:
		ltlen.append(len(i))
	lnm = max(ltlen) + 8     #c'est la longueur du nom le plus long plus 8 espaces
	
	
	
	#gestion des symboles * de similitudes
	symb = ""
	
	for i , j in enumerate(sa):
		
		count = 0
		
		for k in ss:
			
			if j != k[i]: count += 1
			
		if count == 0: symb += "*"
		else: symb += " "
	
	#gestion des tirets
	tir = "----:----|" * int(width/10)
	
	

	d = start
	
	f = d + width
	
	
	#initiation de la str devant contenir la repsrésentation de l'alignement
	txt = """
	
	------------------------------------------------------------------------
	                   Genix alignement - format anagene
		
	Alignement des sequences:
	
		{0}
		{1}
		
	Alignement réalisé avec l'algorithme {2}
	
	Longueur de l'alignement: {3}
	
	% d'identités ou matrice de scores:
	
	{4}
	
	
	Legendes: _ deletion, - similitude avec la référence, * similitude
	-------------------------------------------------------------------------
	
	\n\n\n""".format(nsa, "\n                ".join(nnn[:len(nnn)-1]),program, align.get_alignment_length(), identity)

	

	while d <= stop:
		
		if f > stop:
			
			tir = tir[:len(sa[d:stop])]
		
		
		
		num = " "*lnm
		
		for k in np.arange(10, width + 10, 10):
				
			num = num + """{0}{1}""".format("         " if k == 10 else "        " if  d+k <= 100 else "       " if d+k <= 1000 else "      " , d + k if d+k <= stop else "")
		
		
		if f <= stop:
			
			txt = txt + " "*lnm + symb[d:f] + "\n"
			txt = txt + nsa + " "*(lnm - len(nsa)) + san[d:f] + "\n"
			
			for nsb, sbn in zip(nnn,sss):
				txt = txt + nsb + " "*(lnm - len(nsb)) + sbn[d:f] + "\n"
				
			txt = txt + "\n"
			txt = txt + " "*lnm + tir + "\n" + num + "\n\n\n\n"
		
		else:
			
			txt = txt + " "*lnm + symb[d:stop] + "\n"
			txt = txt + nsa + " "*(lnm - len(nsa)) + san[d:stop] + "\n"
			
			for nsb, sbn in zip(nnn,sss):
				txt = txt + nsb + " "*(lnm - len(nsb)) + sbn[d:stop] + "\n"
				
			txt = txt + "\n"
			txt = txt + " "*lnm + tir + "\n" + num + "\n\n\n\n"
			
			
		
		d, f = d + width, f + width



	return txt
	







	
	
	
	
	
	
	

	



class Emboss_aln_reader:
	
	
	"""
	Permet de lire un alignement emboss et d'en extraire les info ainsi que l'alignement
	"""
	
	def __init__(self, emboss_alignement_path):
		
		self.path = emboss_alignement_path
			
		self.aln = AlignIO.read(self.path, "emboss")
				
		self.info = self.chek_info()
		
		
	@property
	def read_aln(self):
		
		with open(self.path, 'rb') as nd:
			
			x = nd.read()
			x = x.decode()
		
		return x
	
	
	def chek_info(self):
		
		"""lit le fichier emboss généré par un alignement needle ou water pour en extraire les informations d'en-tête
		retourne un dictionnaire d'info"""
	
		
		info = {}
		
	
		with open(self.path,"rb") as emb:
	
			lines = emb.readlines()
	
	
			for i in lines:
				
				i = i.decode()
	
				if i[0] == "#":
	
	
					#recherche du type de l'alignement needle ou water
					if "Program" in i:
	
						l = i.split(' ')
		
						t = l[-1]				
	
						info['program'] = t.strip('\n')
					
					#recherche du nom de la sequence 1
					elif "# 1: " in i:
	
						l = i.split(' ')
						t = l[-1]
						info['name1'] = t.strip('\n')
	
					#recherche du nom de la sequence 2
					elif "# 2: " in i:
				
						l = i.split(' ')
						t = l[-1]
						info['name2'] = t.strip('\n')
	
					#recherche du nom de la matrice de substitution
					elif "# Matrix: " in i:
	
						l = i.split(' ')
						t = l[-1]
						info['matrix'] = t.strip('\n')
	
					#recherche de la valeur de la penalité d'ouverture (gap)
					elif "# Gap_penalty: " in i:
	
						l = i.split(' ')
						t = l[-1]
						info['gap_penalty'] = float(t.strip('\n'))
						
						
					#recherche de la valeur de la penalité  d'extension d'ouverture (gap)
					elif "# Extend_penalty: " in i:
	
						l = i.split(' ')
						t = l[-1]
						info['extend_penalty'] = float(t.strip('\n'))
						
						
					#recherche de la longueur de l'alignement
					elif "# Length: " in i:
	
						l = i.split(' ')
						t = l[-1]
						info['length'] = float(t.strip('\n'))
					
					
					#recherche du pourcentage d'identités
					elif "# Identity: " in i:
	
						l = i.split(' ')
						t = l[-1]
						t = t.strip('\n')
						info['identity'] = float(t[1:-2])
					
					
					#recherche du pourcentage de similarités
					elif "# Similarity: " in i:
	
						l = i.split(' ')
						t = l[-1]
						t = t.strip('\n')
						info['similarity'] = float(t[1:-2])
	
					#recherche du pourcentage de gap
					elif "# Gaps: " in i:
	
						l = i.split(' ')
						t = l[-1]
						t = t.strip('\n')
						info['gaps'] = float(t[1:-2])
						
						
					#recherche du score de l'alignement
					elif "# Score: " in i:
	
						l = i.split(' ')
						t = l[-1]
						info['score'] = float(t.strip('\n'))
						
						
		return info
		
		
	
	@property
	def program(self):
		
		return self.info['program']
				
				
	@property
	def name1(self):
		
		return self.info['name1']
				
	@property
	def name2(self):
		
		return self.info['name2']
		
	@property
	def matrix(self):
		
		return self.info['matrix']
		
	@property
	def gap_penalty(self):
		
		return self.info['gap_penalty']

	@property
	def extend_penalty(self):
		
		return self.info['extend_penalty']
		
	@property
	def length(self):
		
		return self.info['length']
		
	@property
	def identity(self):
		
		return self.info['identity']
		
	@property
	def similarity(self):
		
		return self.info['similarity']
		
	@property
	def gaps(self):
		
		return self.info['gaps']
		
	@property
	def score(self):
		
		return self.info['score']
				
	@property
	def seq1(self):
		
		return self.aln[0].seq
		
	@property
	def seq2(self):
		
		return self.aln[1].seq
		
	@property
	def description1(self):
		
		return self.aln[0].description
		
	@property
	def description2(self):
		
		return self.aln[1].description
		
		
		
		


def showgenix(align, start = 0, stop = None, width = 60, program = None, identity = None):
	
	
	"""présentation d'un alignement nucléotique ou peptidique à une lettre decouppé en tronçons avec regle graduée.
	start et stop doivent etre des entiers; width doit etre un multiple de 10
	
	Arguments:
		align: l'alignement (str) retourné par la fonction AlignIO.read de biopython: exemple : align = AlignIO.read(path, "clustal")
		
		start: (int) index du premier nucléotide représenté
		stop: (int) index du dernier nucléotide représenté
		width: (int) largeur des lignes; width doit etre un multiple de 10
		program: (str) nom du programme ayant réaliser l'alignement (facultatif)
		identity: (str, int ou float) % d'identité ou matrice de score 
	"""
	
	
	s = list(align)
	
	#sequence de reférence 
	sa = str(s[0].seq)
		
	#autres sequence de l'alignement
	ss = [str(i.seq) for i in s[1:]]



	#verification des options

	if width/10 != int(width/10):
		return "Value Error: le parametre width doit etre un multiple de 10"

	if stop == None: 
		stop = len(sa)

	elif type(stop) == type(int()) and stop <= len(sa): 
		pass

	else: 
		return "Value Error: le parametre stop doit etre un entier inférieur ou égal à la longueur de l'alignement"

	if type(start) != type(int()):
		return "Value Error: le parametre start doit etre un entier"


	#récupération des noms des séquences
	
	#nom de la séquence de référence
	nsa = str(s[0].id)
	
	#liste des noms des autres séquences de l'alignement
	nnn = [str(i.id) for i in s[1:]]
	
	
	#recherche du nom de séquence le plus long et gestion des espaces entre nom et séquences
	lt = nnn
	lt.append(nsa)
	
	ltlen = []
	for i in lt:
		ltlen.append(len(i))
	lnm = max(ltlen) + 8     #c'est la longueur du nom le plus long plus 8 espaces
	
	
	
	#gestion des symboles * de similitudes
	#symb = align._star_info ca ne marche pas toujours !!!
	symb = ""
	
	for i , j in enumerate(sa):
		
		count = 0
		
		for k in ss:
			
			if j != k[i]: count += 1
			
		if count == 0: symb += "*"
		else: symb += " "
	
	
	#gestion des tirets
	tir = "----:----|" * int(width/10)
	
	
	#initiation des intervals
	d = start
	
	f = d + width
	
	
	
	#initiation de la str devant contenir la repsrésentation de l'alignement
	txt = """
	
	------------------------------------------------------------------------
	                    Genix alignement - format genix
		
	Alignement des sequences:
	
		{0}
		{1}
		
	Alignement réalisé avec l'algorithme {2}
	
	Longueur de l'alignement: {3}
	
	% d'identités ou matrice de scores:
	
	{4}
	
	Legendes: - deletion, * similitude
	-------------------------------------------------------------------------
	
	\n\n\n""".format(nsa, "\n                ".join(nnn[:len(nnn)-1]),program, align.get_alignment_length(), identity)

	

	while d <= stop:
		
		if f > stop:
			
			tir = tir[:len(sa[d:stop])]
		
		
		
		num = " "*lnm
		
		for k in np.arange(10, width + 10, 10):
				
			num = num + """{0}{1}""".format("         " if k == 10 else "        " if  d+k <= 100 else "       " if d+k <= 1000 else "      ", d + k if d+k <= stop else "")
		
		
		if f <= stop:
			
			txt = txt + " "*lnm + symb[d:f] + "\n"
			txt = txt + nsa + " "*(lnm - len(nsa)) + sa[d:f] + "\n"
			
			for nsb, sb in zip(nnn,ss):
				txt = txt + nsb + " "*(lnm - len(nsb)) + sb[d:f] + "\n"
				
			txt = txt + "\n"
			txt = txt + " "*lnm + tir + "\n" + num + "\n\n\n\n"
		
		else:
			
			txt = txt + " "*lnm + symb[d:stop] + "\n"
			txt = txt + nsa + " "*(lnm - len(nsa)) + sa[d:stop] + "\n"
			
			for nsb, sb in zip(nnn,ss):
				txt = txt + nsb + " "*(lnm - len(nsb)) + sb[d:stop] + "\n"
				
			txt = txt + "\n"
			txt = txt + " "*lnm + tir + "\n" + num + "\n\n\n\n"
			
			
		
		d, f = d + width, f + width



	return txt
	
	




def treeliste(s):
	
	"""convertion d'une str representant une seq petpitidique à 3 lettres : MetValGlu en une liste d'acides aminés ['Met','Val','Glu']"""
	
	
	d = 0
	f = 3
	
	tl = []
	

	while f <= len(s):
		
		tl.append(s[d:f])
		d , f = f , f+3

	return tl




def showgenixp3(align, start = 0, stop = None, width = 30, program = None, identity = None):
	
	
	"""présentation d'un alignement peptidique avec langage 3 lettres decouppé en tronçons avec regle graduée.
	start et stop doivent etre des entiers; width doit etre un multiple de 10
	
	Arguments:
		align: l'alignement (str) retourné par la fonction AlignIO.read de biopython: exemple : align = AlignIO.read(path, "clustal")
		
		start: (int) index du premier nucléotide représenté
		stop: (int) index du dernier nucléotide représenté
		width: (int) largeur des lignes; width doit etre un multiple de 10
		program: (str) nom du programme ayant réaliser l'alignement (facultatif)
		identity: (str, int ou float) % d'identité ou matrice de score 
	"""
	
	
	s = list(align)
	
	#sequence de reférence 
	sa = str(s[0].seq)
		
	#autres sequence de l'alignement
	ss = [str(i.seq) for i in s[1:]]
	
	
	
	#verification des options

	if width/10 != int(width/10):
		return "Value Error: le parametre width doit etre un multiple de 10"

	if stop == None: 
		stop = len(sa)

	elif type(stop) == type(int()) and stop <= len(sa): 
		pass

	else: 
		return "Value Error: le parametre stop doit etre un entier inférieur ou égal à la longueur de l'alignement"

	if type(start) != type(int()):
		return "Value Error: le parametre start doit etre un entier"
	
		
	#passage au langage a 3 lettres
	sa = seq3(sa, custom_map={"*": ""}, undef_code=' - ')
	
	#transformation en list d'acides amines ['Met','Val','Glu']
	sa = treeliste(sa)
	
	#les deletions notées - dans l'alignement sont converties par seq3 en Xaa , il faut donc les corriger.... inutile avec undef_code=' - '
	#for i, j in enumerate(sa): 
		
		#if j == 'Xaa':
			
			#sa[i] = ' - '
	
	
	
	#creation d'une liste de liste des autres sequences avec langage a 3 lettres
	sss = []
	
	for i in ss:
		
		j = seq3(i, custom_map={"*": ""}, undef_code=' - ')
		
		ssl = treeliste(j)
		
		#for i, j in enumerate(ssl): 
		
			#if j == 'Xaa':
			
				#ssl[i] = ' - '
		
		sss.append(ssl)
		
	ss = sss
	
		
	#gestion des tirets
	tir = " -  -  -  -  :  -  -  -  -  | " * int(width/10)
	
	
	#récupération des noms des séquences
	
	#nom de la séquence de référence
	nsa = str(s[0].id)
	
	#liste des noms des autres séquences de l'alignement
	nnn = [str(i.id) for i in s[1:]]
	
	
	#recherche du nom de séquence le plus long et gestion des espaces entre nom et séquences
	lt = nnn
	lt.append(nsa)
	
	ltlen = []
	for i in lt:
		ltlen.append(len(i))
	lnm = max(ltlen) + 8     #c'est la longueur du nom le plus long plus 8 espaces
	
	
	
	#gestion des symboles * de similitudes
	symb = ""
	
	for i , j in enumerate(sa):
		
		count = 0
		
		for k in ss:
			
			if j != k[i]: 
				count += 1
			
		if count == 0: 
			symb += " * "
			
		else: 
			symb += "   "
	
	symb = treeliste(symb)
	
	
	
	#initiation des intervals
	d = start
	
	f = d + width
	

	#initiation de la str devant contenir la repsrésentation de l'alignement
	txt = """
	
	------------------------------------------------------------------------
	                    Genix alignement - format genix
		
	Alignement des sequences:
	
		{0}
		{1}
		
	Alignement réalisé avec l'algorithme {2}
	
	Longueur de l'alignement: {3}
	
	% d'identités ou matrice de scores:
	
	{4}
	
	Legendes: - deletion, * similitude
	-------------------------------------------------------------------------
	
	\n\n\n""".format(nsa, "\n                ".join(nnn[:len(nnn)-1]),program, align.get_alignment_length(), identity)

	
	while d <= stop:
		
		if f > stop:
			
			tir = tir[:len(''.join(sa[d:stop]))]
		
		
		num = " "*lnm
		
		for k in np.arange(10, width + 10, 10):
				
			num = num + """{0}{1}""".format("                            " if k == 10 else "                            " if  d+k <= 100 else "                           " , d + k if d+k <= stop else "")
		
		
		if f <= stop:
			
			txt = txt + " "*lnm + "".join(symb[d:f]) + "\n"
			txt = txt + nsa + " "*(lnm - len(nsa)) + "".join(sa[d:f]) + "\n"
			
			for nsb, sb in zip(nnn,ss):
				txt = txt + nsb + " "*(lnm - len(nsb)) + "".join(sb[d:f]) + "\n"
				
			txt = txt + "\n"
			txt = txt + " "*lnm + tir + "\n" + num + "\n\n\n\n"
		
		else:
			
			txt = txt + " "*lnm + "".join(symb[d:stop]) + "\n"
			txt = txt + nsa + " "*(lnm - len(nsa)) + "".join(sa[d:stop]) + "\n"
			
			for nsb, sb in zip(nnn,ss):
				txt = txt + nsb + " "*(lnm - len(nsb)) + "".join(sb[d:stop]) + "\n"
				
			txt = txt + "\n"
			txt = txt + " "*lnm + tir + "\n" + num + "\n\n\n\n"
			
					
		d, f = d + width, f + width

	return txt
	
	
def supraccent(x):
	
	y = unicodedata.normalize('NFKD', x).encode('ascii', 'ignore')
	
	return y.decode()
	
def edi2fasta(path_edi, path_fas):
	
	"""
	Transforme un fichier edi en fasta
	
	Argument:
		
		path_edi: path du fichier en .edi
		path_fasta: path du fichier en .fasta
		
		les extensions sont obligatoires !!!
		
	"""

	#verification des extensions
	if path_edi[len(path_edi)-4:] != '.edi':
		return 'Error path_edi extension different de .edi'
		
		
	if path_fas[len(path_fas)-4:] != '.fas':
		return 'Error path_fas extension different de .fas'
		
		
	#ouverture et lecture du fichier .edi
	try:
		with open(path_edi,'r', encoding='ISO-8859-1') as f:
			r = f.read()
			r = supraccent(r)
			
			
	except FileNotFoundError:
		return 'fichier edi inexistant'

	#decoupage en liste de sequences (si plusieurs sequences dans le edi)
	l = r.split(';-')
	
	
	list_fasta = []
	
	for i in l[:len(l)-1]:
		
		#decoupage de la sequence
		a = i.split('\n;')
		
		try:
			name = a[1][1:] #recuperation du nom de la sequence et suppression de l'espace
			if '-' in name:
				name = name.replace('-', '_')
			
			
			if len(a) == 6:
				des = a[4][1:] #description de la sequence
			
			elif len(a) > 6:
				des = a[4:len(a)-1]
				des = "".join(des)
			else:
				des = ''
			
			seq = a[len(a)-1] #la seq est en dernier dans la liste
			name = name.replace(" ","_")
			
			
		except IndexError:
			return 'Le fichier .edi presente une structure anormale'
		

		else:
			#suppression des espaces et des \n de la sequence
			y = list(seq) #conversion str to list
			while ' ' in y:
				y.remove(' ')
			while '\n' in y:
				y.remove('\n')
				
			seq = "".join(y) #reconversion list to str
			
			fas = SeqRecord(Seq(seq), id = name , description=des)
			
			list_fasta.append(fas)
			
	
	#ecriture du fasta
	SeqIO.write(list_fasta, path_fas, "fasta")
	
		
