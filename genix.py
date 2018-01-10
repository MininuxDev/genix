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
from genixm import licence

import argparse
import textwrap









if os.getenv('LANG') == "fr_FR.UTF-8":	
	
	parser = argparse.ArgumentParser(prog = 'Genix', description= "Outils d'analyses génétiques et phylogénétiques")
	
else:
	
	parser = argparse.ArgumentParser(prog = 'Genix', description= "Genetic tools for education")
	
	
	
parser.add_argument('--version', action='version', version='%(prog)s 0.1')

parser.add_argument('--license', action='store_true')

parser.add_argument('--edi2fas', nargs = 2, help="transforme un edi en fasta le path du fichier edi doit se terminer en .edi tandis que le path du fichier fasta doit se terminer en .fas" )



if os.getenv('LANG') == "fr_FR.UTF-8":

	subparsers = parser.add_subparsers(help="genix - outils d'analyses genetiques pour l'éducation", dest = 'sub')

else:
	
	subparsers = parser.add_subparsers(help="genix - genetic tools for education", dest = 'sub')


#menu

parser_menu = subparsers.add_parser('menu')




args = parser.parse_args()



print('\n\n')





	

if args.license:
	
	t = textwrap.dedent("""\

	    =====================================================================

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
	
	    =====================================================================

	""")
	
	print(t)
	
	
if args.edi2fas:
	
	p = args.edi2fas
	
	edi = p[0]
	
	fasta = p[1]
		
	edi2fasta(edi,fasta)
	
		
				

if args.sub == "menu": 
	
	licence()
