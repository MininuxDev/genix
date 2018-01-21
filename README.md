# Presentation de Genix.

Genix est un programme de génie génétique principalement destiné à l’étude des SVT. Ce programme s’utilise dans un terminal GNU/Linux et comporte un menu interactif permettant des réaliser l’ensemble des traitements possibles. Il a été développé dans le but d’offrir une alternative aux logiciels (libres et non libres) pédagogiques de traitement de séquences déjà existant mais utilisables uniquement sous Windows.

## Fonctionnalités.

- Recherche et ajout de séquences (nucléotidiques et peptidiques)
- Création de séquences (nucléotidiques et peptidiques)
- Transcritpion
- Traduction
- Rétro-transcription
- Affichage d’informations sur les séquences
- Alignements et comparaisons simples locales et globales 
- Alignements et comparaisons multiples (Clustalw)
- Choix du code génétique utlilisé pour la traduction
- Représentation des tables des différents codes génétiques
- Affichage des séquences peptidiques par une représentation à une ou trois lettres par acide aminé
- Choix de la représentations des alignements (genix, emboss, clustalw, phylip, anagène-like)
- Sauvegarde des séquences créées, ou modifiées par transcription , traduction ou rétro-transcription
- Conversion de séquences au format .edi en fasta
- Création de fasta mono et multi-séquences
- Affichage des scores et/ou matrices de similitudes associés aux alignements

## Installation.

Genix est un programme qui fonctionne uniquement sous GNU Linux. Il est écrit en python et requiert une version de python >= 3.4

### Installation des dépendances.

Genix requiert des dépendances python qui ne sont pas présentes dans la librairie standard (numpy, biopython, xhtml2pdf, pandas) et des dépendances non python (Clustalw, Emboss, Rebase, Phylip).

#### Installation des dépendances Python.

	pip install numpy biopython xhtml2pdf pandas


#### Installation de Clustalw

http://www.clustal.org/clustal2/

Sur Debian:

	sudo apt-get install clustalw

#### Installation de Emboss

http://emboss.sourceforge.net/download/

Sur Debian:

	sudo apt-get install emboss

#### Installation rebase pour emboss. 

Rebase est indispensable pour utiliser les fonctions liées aux enzymes de restriction; il s'agit notamment des fonctions restrict et remap.

1: placez vous dans le répertoire de rebase

	cd /usr/share/EMBOSS/data/REBASE

2: télécharger les fichiers withrefm et proto depuis le serveur ftp de rebase

	sudo wget ftp://ftp.ebi.ac.uk/pub/databases/rebase/withrefm*

Puis 

	sudo wget ftp://ftp.ebi.ac.uk/pub/databases/rebase/proto*

3: décompresser 

	sudo uncompress withrefm.707.gz

Et 

	sudo uncompress proto.707.gz

Vous n'aurez pas forcément le version 707.... c'est à adapter

4: extraire rebase

	rebaseextract

Et voilà... c'est terminé.

#### Installation de Phylip

http://evolution.genetics.washington.edu/phylip/install.html


Sur Debian:

	sudo apt-get install phylip

### Installation de genix

Téléchargement avec wget:

	~ $ wget https://github.com/YannBouyeron/genix/archive/master.zip
	
	~ $ unzip master.zip
	
	~ $ mv genix-master genix

Telechargement avec git:

	~ $ git clone https://github.com/YannBouyeron/genix.git

Rendre le programme executable:

	~ $ cd genix && mv genix.py genix && sudo chmod 755 genix && sudo cp genix* /usr/local/bin && cd ..

## Utilisation

Accéder au menu interactif:

	~ $ genix -m

Accéder au menu interactif (si vous n’avez pas rendu le programme exécutable):

	~ $ cd genix
	~/genix $ python3 genix.py -m

Accéder à l’aide:

	~ $ genix -h

## Désinstallation

	~ $ cd /usr/local/bin/
	/usr/local/bin $ sudo rm genix*

