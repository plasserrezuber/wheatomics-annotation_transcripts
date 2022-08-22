#!/bin/bash

#connexion 7srv pui topaze
#verif droits sur dossier courant (.) (parent etant le ..) en faisant ls -la : j'ai bien les droits en ecriture
#verif taille des directory a transferer de topaze vers hpc2 avec "df -sch ./*"
#verif espace dispo sur 7srv avec "df -h"

#tentative pour eviter de taper le password pour chacun des 80 directory:
#export RSYN_PASSWORD="password"
#rsync -auv rsync://lasserrp@topaze.ccc.cea.fr:/ccc/store/cont007/fg0118/fg0118/rawdata/projet_CGC/${line} /home/plasserre/data/wheatomics/rnaseq/
#echec connexion a topaze avec le rsync:// protocole

###autre idee: faire un dossier sur topaze contenant des liens symbolique pointant vers les dossiers a transferer.

#transfert liste des dossier a transferer depuis 7srv
scp /home/plasserre/data/wheatomics/rnaseq/list_wheatomics_rnaseq_grain_5cv.txt lasserrp@topaze.ccc.cea.fr:/ccc/store/cont007/fg0118/fg0118/rawdata/projet_CGC/

#sur topaze:
mkdir rnaseq_grain_5cv
mv list_wheatomics_rnaseq_grain_5cv.txt /ccc/store/cont007/fg0118/fg0118/rawdata/projet_CGC/rnaseq_grain_5cv/

while read line;
do
    ln -s /ccc/store/cont007/fg0118/fg0118/rawdata/projet_CGC/${line} /ccc/store/cont007/fg0118/fg0118/rawdata/projet_CGC/rnaseq_grain_5cv/${line}
done < /ccc/store/cont007/fg0118/fg0118/rawdata/projet_CGC/rnaseq_grain_5cv/list_wheatomics_rnaseq_grain_5cv.txt


#sur 7srv:
#screen
screen -S topTO7

#puis:
#-L, --copy-links            transform symlink into referent file/dir
rsync -rL lasserrp@topaze.ccc.cea.fr:/ccc/store/cont007/fg0118/fg0118/rawdata/projet_CGC/rnaseq_grain_5cv /home/plasserre/data/wheatomics/rnaseq/

#pour detacher le screen sans fermer le terminal virtuel
#Ctrl+a d

#pour lister
screen -ls
#pou rvoir la progression du transfert:
#-c: taille totale
du -ch ~/data/wheatomics/rnaseq/rnaseq_grain_5cv

#pour rattacher le screen 
screen -r topTO7

##puis de 7srv vers hpc2
rsync -rv /home/plasserre/data/wheatomics/rnaseq/rnaseq_grain_5cv hpc2:/storage/groups/gdec/shared/triticum_aestivum/wheatomics/rnaseq/

##kill screen: exit