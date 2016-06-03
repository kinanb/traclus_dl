# traclus_dl

TraClus_DL est un programme qui permet d'identifier les corridors à partires des données origine-destination (lignes de désir), TraClus-DL inspiré du programme Traclus (Clustering Trajectoire, publié par Lee et al en 2007).
Les données d'entrées comprennent : 
1)	Identification unique pour chaque ligne de désir (ID), 
2)	Poids de la ligne de désir (FEX) ou un facteur de pondération. En cas d’absence d’un facteur de pondération la valeur (1) doit être utilisée ;
3)	Coordonnée x d’origine (x_o) ; 
4)	Coordonnée y d’origine (x_o) ; 
5)	Coordonnée x de destination (x-d) ;
6)	 Coordonnée y de destination (y_d).
Le fichier contenant les données d’entrées (filename) est un fichier texte. Chaque ligne du fichier décrit une ligne de désir : 
<ID> <FEX> <x_o> <y_o> <x_d> <y_d>
Les données de sortie sont récapitulées dans deux fichiers : 
  1) un fichier contenant l'ensemble des corridors identifiés. Ces corridors peuvent être facilement visualisés dans un logiciel tel que QGIS) ;
  2) un fichier reliant chaque corridor aux lignes de désir associées.
Le logiciel est écrit en python. Il est exécutable à partir de la ligne de commande où Python est installé comme :
python Traclus_DL.py <filename> <max_dist> <min_density> <angle_max> <segment_size>

Plus de détails sur les paramètres d'entrées et le programme sont disponible dans : 
Bahbouh, K., Wagner, J. R., Morency, C., & Berdier, C. (2015). TraClus-DL: A Desire Line Clustering Framework to Identify Demand Corridors. In Transportation Research Board 94th Annual Meeting (No. 15-3508).
https://trid.trb.org/view.aspx?id=1338139
