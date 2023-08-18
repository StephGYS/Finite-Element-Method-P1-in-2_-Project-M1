
--> proj.ffem
Maillage non structuré d'un domaine en forme de "fusée". 
Nous définissons la distance: R( le rayon de la coiffe) et les sommets caractéristiques avec "real x=xi, y=yi; avec (xi,yi) 
la coordonnée du point" afin de définir la forme de "fusée".
Ensuite on paramétrise les bords avec la fonction "border name(t=deb, fin){x=x(t) ;y=y(t) ; label =num_label} avec x(t)=(1-t).x1+(t).x2 et y(t)=(1-t).y1+t.y2" ,
puis on discrétise le domaine avec la fonction "mesh mon_maillage=buildmesh(ce(n)+a(n)+b(n)+c(n)+d(n)+e(n)+f(n)+g(n)); avec ce: la coiffe, a,b,c,d,e,f,g: les segments de la fusée et n: le nombre de nœuds. 
Finalement pour l'affichage, on se sert de "plot(nom_maillage)" et de "savemesh(nom_maillage,"nom_fichier.msh");" pour une éventuelle sauvegarde.    


-->validation.ffem
-PRETRAITEMENT: Après avoir lire le fichier ".msh", nous avons définit les parametres(ka,Kb,e,V,x0,y0,alphaD,uE,uD....),la fonction "f" et la condition de Fourier/Robin (ici uf=uE) 
   -TRAITEMENT: Nous passons à la résolution. D'abord on définit le type d'espace d'élément fini "Vh" avec "fespace", les inconnues de Vh "uh et vh". Ensuite notre équation à résoudre "a(uh,vh)=l(vh)" avec "PROBLEM(.....,solver=Cholesky)", on ajoute ensuite la condition de Dirichlet avec "on(RefDir, uh=uD)" et la condition de Fourier avec "on(RefFou, uh=uF)"
   -Affichage de la solution: nous définissons le type de forme variationnelle "b(uh,vh)" grâce à "varf" que nous sauvegardons dans la matrice "K". 
 Remarques : Pour les intégrales simples "int1d()" et doubles "int2d".
 


-->projj.py
   1)"""Pre-traitement"""  
Il s'agit dans cette partie de définir  les informations  nécessaires  à  la résolution  de notre problème  de formulation variationnelle. D’abord  nous avons besoin de charger le maillage et de connaître  son pas et la qualité  afin de l’évaluer : 
--charge_et_affiche_maillage(): Cette fonction permet de lire et afficher un fichier maillage .msh, c'est l'ensemble de la fonction lit_fichier_msh et des 2 procédures définies dans le tp1. 
--pas_et_qualite_maillage(): Cette fonction permet de calculer le pas et la qualité du maillage à partir des coordonnées et du nombre de coordonnées. Elle utilise le théorème de Pythagore pour calculer le pas (max des longueurs des triangles) qui n'est autre que l'hypoténuse. Elle prend la distance du plus petit découpage sur l'axe des abscisses et la plus petite coordonnées de l'axe des ordonnées. 
La qualité est obtenue à partir de la formule vue en cours et d'une deuxième application du théorème de Pythagore pour le calcul de la hauteur du triangle dans le calcul de l'aire du triangle. Ensuite, nous définissons les fonctions  nécessaires  à  la résolution du problème : 
--fct_uE(x,y) : la température extérieure uE 
--fct_f(x,y) : la fonction source donnée  dans le tp2 
--fct_kappa(x,y) :la fonction  de conductivité notée k 
--fct_alpha(x,y) : le facteur de transfert noté alpha 
-- fct_beta:le facteur de transtert sur dirichlet
--coeffelem_p1_rigid(x1,x2,x3,y1,y2,y3,mes,nbe,coord,tri) la matrice de rigidité de l’élément Tl notée k. 
les paramètres  d’entrée sont : Les coordonnées des nœuds, l’aire de triangle,  le nombre de nœuds,  les coordonnées des nœuds et la table des indices de chaque triangle(maillage Tl). Puis, en utilisant les formules données, nous retournons k. 
--masse(x1,x2,x3,y1,y2,y3,mes,nbe,coord,tri) la matrice de masse notée ml. la matrice de rigidité notée k. 
les paramètres  d’entrée sont : Les coordonnées des nœuds, l’aire de triangle,  le nombre de nœuds,  les coordonnées des nœuds et la table des indices de chaque triangle(maillage Tl). Puis, en utilisant la formule donnée, nous retournons ml. 
--coeffelem_p1_source(x1,x2,x3,y1,y2,y3,mes,nbe,coord,tri) le vecteur source de l’élément Tl noté fl.
 les paramètres  d’entrée sont : Les coordonnées des nœuds, l’aire de triangle,  le nombre de nœuds,  les coordonnées des nœuds et la table des indices de chaque triangle(maillage Tl). 
Puis, en utilisant la formule donnée, nous retournons fl. 
--coeffelem_P1_transf(x1,x2,y1,y2,mes,nba,coord,ar) Le vecteur pour l’arête Aa noté ea. les paramètres  d’entrée sont : Les coordonnées des nœuds, l’aire de triangle,  le nombre de nœuds,  les coordonnées des nœuds et la table des indices des sommets  de chaque arêtes. Puis, en utilisant la formule donnée, nous retournons ea. 
--coeffelem_P1_poids(x1,x2,y1,y2,mes,nba,coord,ar) La matrice pour l’arête Aa notée pa. les paramètres  d’entrée sont : Les coordonnées des nœuds, l’aire de triangle,  le nombre de nœuds,  les coordonnées des nœuds et la table des indices des sommets  de chaque arêtes. Puis, en utilisant la formule donnée, nous retournons pa.  
--coeffelem_P1_transf_dir(x1,x2,y1,y2,mes,nba,coord,ar)
--coeffelem_P1_poids_dir(x1,x2,y1,y2,mes,nba,coord,ar)
   2)"""traitement"""
 Le but de cette partie est de trouver la solution  approchée Uh. N’oublions pas que le système  à  résoudre  est AU=F, nous procéderons d’abord à l’assemblage de la MEF pour avoir A,k(matrice de rigidité)et F. 
--assemblage_EF_P1(nbn,nbe,nba,tri,ar,refa) ; les paramètres  d’entrée sont : le nombre de nœuds, le nombre d’éléments(maillage Tl),le nombre total d’arêtes de frontière, la table des indices des sommets  de chaque éléments, la table des indices des sommets  de chaque arêtes et les indices des références  pour identifier  les arêtes de bord. Puis, en utilisant la script donné, nous retournons A,k et F. Il faut remarquer que nous avons utilisé la formule de l’aire d’u triangle vu au TD1, aussi vu la condition de Fourier/Robin (alpha grand=> Dirichlet ), les termes non nuls de A sont ceux où  les sommets  des arêtes de bord appartiennent au bord du domaine(i.e refa[a,0]==1). 
--facto_chol(A) solveur cholesky AU=F, pourquoi ? Vu que A est sdp,il est idéal  de résoudre le système  avec la méthode de cholesky. Dans cette fonction , on décompose A tel que A=CC^t avec C ne matrice  triangulaire  inférieure. On sauvegarde donc C et C^t. 
--descente(C,F) Le système  revient  donc à (CC^t)U=F alors dans cette fonction , nous  résolvons Cy=F en descendant pour avoir y qu’on sauvegarde. 
--remontee(U,y) Cela reviendra  finalement à  résoudre (C^t)U=y en remontant pour avoir U. On sauvegarde  finalement  U dans la variable Uh (solution  approchée).  

  3)"""Post-traitement""" 
Le but de cette partie est de pourvoir valider et interpréter nos résultats  obtenus, nous vérifions en premier la convergence de Uh : 
--conv(nbn) Pour vérifier  la convergence, nous traçons les fonctions chapeaux grâce à l’instruction plot_trisurf(coord[ :,0],coord[ :,1],tri,Uh).N’oublions pas que les composantes de Uh dans la base des fonctions chapeaux  sont les valeurs de Uh prises à chaque nœuds.  Cette fonction renvoie donc le graphe des fonctions  chapeaux  et les composantes de la solution exacte U. En second, nous calculer l’erreur en semi-norme: 
--eh, en semi-norme H^1 En se servant de la formule de l’erreur donnée dans le tp, nous avons finalement  l’erreur vu que nous avons déjà Uh, U et k(matrice  de rigidité).
 En dernier : nous traçons  la forme de la matrice A grâce  à  l’instructions spy(A), le graphe de la solution approchée Uh. 

  4)"""test Affichage des résultats numériques""" Cette partie sert de test pour les précédentes  parties, nous avons exécuté les instructions  précédentes  en chargeant le fichier maillage "fusee.msh".   
