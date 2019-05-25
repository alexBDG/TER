---------------------------------------
Utilisation du program avec "data.txt"
---------------------------------------

CONDITION D'UTILISATION
  - Plaque rectangulaire
  - Numérotation des cotés dans le *.mesh :
      _ Ref 1 : Bord bas
      _ Ref 2 : Bord droit
      _ Ref 3 : Bord haut
      _ Ref 4 : Bord gauche

---------------------------------------
---------------------------------------

DETAILS DES PARAMETRES MODIFIABLES
  - Choix du maillage

  - Temps initial/final et pas de Temps

  - Coefficients propres à l'equation de propagation

  - Frequence a laquelle les solutions sont enregistrées

  - Conditions intiniales
      _ Center : Déplacement initial définit sur un cercle
      de centre (x0, y0) et de rayon r.
      _ Rectangular : Deplacement positif selon x défini sur
      une bande verticale de coordonnée y0 et d'eppaisseur thn
      _ no : Conditions initiales nulles

  - Conditions aux bords
      _ Neumann : Aucune modification est faite dans le programme
      _ Dirichlet_h : Dirichlet homogène de valeurs A et d'orientation
      dépendant du choix [ c (cisaillement) / tc (traction compression) ]
      _ Dirichlet_nh : Dirichlet non homogène de la forme
      A.cos(w*2*pi*t)*exp((t-tm)^2)

  - Excitation sur une partie ou sur tout le bord gauche et/ou haut
  de la forme A.cos(w*2*pi*t)*exp((t-tm)^2) et d'orientation dépendant
  du choix [ c (cisaillement) / tc (traction) ]
  REMARQUE : A activer seulement si il n'y a pas déjà des conditions
  de Dirichlet non homogènes sur le meme bord.

  - Possibilitée d'activer au plus 2 capteurs placés en (x0,y0) et
  de taille précisable au niveau de "precision"

  - Non du fichier dans lequel sont enregistrées les solutions
