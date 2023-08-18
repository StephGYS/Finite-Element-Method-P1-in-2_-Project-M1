from pylab import *
from mpl_toolkits.mplot3d.axes3d import Axes3D
"""Pre-traitement"""
# definition de la fonction qui affiche les coord d'un maillage
def charge_et_affiche_maillage(f):
    # lire une ligne 
    ligne = f.readline() #// ou faire un f.(touche tab) // print(ligne)
    data = ligne.split() # pour avoir une liste de chaine de charactère
    nbn = int(data[0])
    nbe = int(data[1])
    nba = int(data[2])
    coord = zeros((nbn,2), double)
    tri = zeros((nbe,3), int)
    ar = zeros((nba,2), int)
    refn = zeros((nbn,1), double)
    reft = zeros((nbe,1), double)
    refa = zeros((nba,1), int)
    
    for i in range(len(coord)):
        # pour lire la ligne suivante on réecrit 
        ligne = f.readline()
        data = ligne.split()
        coord[i,0] = double(data[0]) # pour lire un reel
        coord[i,1] = double(data[1]) # pour lire un reel
        refn[i] = double(data[2])
    
    for i in range(len(tri)):
        ligne = f.readline()
        data = ligne.split()
        tri[i,0] = double(data[0]) # pour lire un reel
        tri[i,1] = double(data[1]) # pour lire un reel
        tri[i,2] = double(data[2])
        reft[i] = double(data[3])
    
    for i in range(len(ar)):
        ligne = f.readline()
        data = ligne.split()
        ar[i,0] = double(data[0]) # pour lire un reel
        ar[i,1] = double(data[1]) # pour lire un reel
        refa[i] = double(data[2])
        
    tri = tri-1
    ar = ar-1
    figure(3)
    triplot(coord[:,0],coord[:,1],tri)
    for i in range(nbn):
        text(coord[i,0],coord[i,1],str(i), bbox=dict(facecolor='red', alpha=0.5))
    figure(4)
    triplot(coord[:,0],coord[:,1],tri)
    for i in range(nbn):
        text(coord[i,0],coord[i,1],str(refn[i,0]), bbox=dict(facecolor='red', alpha=0.5))
     
    return nbn, nbe,nba, coord,tri,ar,refn,reft,refa

## definition de la fonction qui calcul le pas et la qualité d'un maillage
def pas_et_qualite_maillage(nbe,coord,tri):
    h=[]
    Qt=[]
    for i in range(nbe):
        x1=coord[tri[i,0],0]
        x2=coord[tri[i,1],0]
        x3=coord[tri[i,2],0]
        y1=coord[tri[i,0],1]
        y2=coord[tri[i,1],1]
        y3=coord[tri[i,2],1]
        L1=sqrt((x1-x2)**2+(y1-y2)**2)
        L2=sqrt((x2-x3)**2+(y2-y3)**2)
        L3=sqrt((x3-x1)**2+(y3-y1)**2)
        L=max(L1,L2,L3)
        h.append(L)
        mes=0.5*abs(((x2-x1)*(y3-y1)-(x3-x1)*(y2-y1)))
        rt=2*mes/(L1+L2+L3)
        qt=sqrt(3)*L/(6*rt)
        Qt.append(qt)
    hmax=max(h)
    q=max(Qt)
    return hmax, q

#la solution exacte u
#def fct_u(x,y):
#    u=1+sin(pi*x/2)+x*(x-4)*cos(pi*y/2)
#    return u

#la temperature exterieure uE
V=500000;  #a faire varier
def fct_uE(x,y):
    uE=(10**(-3))*(V**2)
    return uE

#la temperature uD sur dirichlet  #donnée
def fct_uD(x,y):
    uD=293.15
    return uD

#la fonction source
x0=0.015;
y0=0.015;
def fct_f(x,y):
    f= -0.5 * exp(-((x-x0)**2+(y-y0)**2) )  #donnée    
    return f

#la fonction de conductivité
def fct_kappa(x,y):
    ka=0.025     #W/mK pour T=293 Kelvin  #donnée 
    return ka

# le facteur de transtert sur fourier
Kb=237; #237W/mK pour T=293 Kelvin   #fixe
e=100;       #a faire varier
def fct_alpha(x,y):
    alpha=Kb/e # au bord de Fourier-Robin
    return alpha

# le facteur de transtert sur dirichlet
def fct_beta(x,y):
    beta=10**(8) # au bord de dirichlet
    return beta

#la matrice de Tl
def coeffelem_p1_rigid(x1,x2,x3,y1,y2,y3,mes,nbe,coord,tri):
    k=zeros((3, 3))
    ka=fct_kappa(1,1)
    # les coordonnées de la matrice
    k[0,0]=ka*((x2-x3)**2+(y2-y3)**2)/(4*mes)
    k[1,1]=ka*((x3-x1)**2+(y3-y1)**2)/(4*mes)
    k[2,2]=ka*((x1-x2)**2+(y1-y2)**2)/(4*mes)
    k[0,1]=ka*(-(x1-x3)*(x2-x3)-(y1-y3)*(y2-y3))/(4*mes)
    k[0,2]=ka*(-(x3-x2)*(x1-x2)-(y3-y2)*(y1-y2))/(4*mes)
    k[1,2]=ka*(-(x2-x1)*(x3-x1)-(y2-y1)*(y3-y1))/(4*mes)
    k[1,0]=k[0,1]
    k[2,0]=k[0,2]
    k[2,1]=k[1,2]
    return k

#la matrice masse
def masse(x1,x2,x3,y1,y2,y3,mes,nbe,coord,tri):
    ml=zeros((3,3))
    #la mesure du triangle
    #les coords du vecteur
    ml[0,0]=mes/3
    ml[1,1]=ml[0,0]
    ml[2,2]=ml[0,0]
    return ml
                                         
#le vecteur de Tl
def coeffelem_p1_source(x1,x2,x3,y1,y2,y3,mes,nbe,coord,tri):
    fl=zeros(3)
    #le milieu du triangle
    mi1=(x1+x2+x3)/3
    mi2=(y1+y2+y3)/3
    fmi=fct_f(mi1,mi2)
    #les coords du vecteur
    for i in range(3):
        fl[i]=mes*fmi/3
    return fl

#le vecteur ea pour l'arret Aa
def coeffelem_P1_transf(x1,x2,y1,y2,mes,nba,coord,ar):
    ea=zeros(2)
    a=fct_alpha(x1,y1)
    #le milieu du triangle
    mi1=(x1+x2)/2
    mi2=(y1+y2)/2
    fmi=fct_uE(mi1,mi2)
    #les coords du vecteur
    for i in range(2):
        ea[i]=mes*a*fmi/2
    return ea

#le vecteur ea_dir pour dirichlet
def coeffelem_P1_transf_dir(x1,x2,y1,y2,mes,nba,coord,ar):
    ea_dir=zeros(2)
    b=fct_beta(x1,y1)
    #le milieu du triangle
    mi1=(x1+x2)/2
    mi2=(y1+y2)/2
    fmi=fct_uD(mi1,mi2)
    #les coords du vecteur
    for i in range(2):
        ea_dir[i]=mes*b*fmi/2
    return ea_dir

# la matrice pa
def coeffelem_P1_poids(x1,x2,y1,y2,mes,nba,coord,ar):
    pa=zeros((2,2))
    a=fct_alpha(x1,y1)
    #les coords du vecteur
    pa[0,0]=mes*a/3
    pa[1,1]=pa[0,0]
    pa[0,1]=mes*a/6
    pa[1,0]=pa[0,1]
    return pa;

#la matrice pa_dir sur dirichlet
def coeffelem_P1_poids_dir(x1,x2,y1,y2,mes,nba,coord,ar):
    pa_dir=zeros((2,2))
    b=fct_beta(x1,y1)
    #les coords du vecteur
    pa_dir[0,0]=mes*b/3
    pa_dir[1,1]=pa_dir[0,0]
    pa_dir[0,1]=mes*b/6
    pa_dir[1,0]=pa_dir[0,1]
    return pa_dir;

"""Traitement"""
#affectation de lamatrice EF
def assemblage_EF_P1(nbn,nbe,nba,tri,ar,refa):
    A=zeros((nbn,nbn))
    F=zeros(nbn)
    k=zeros((nbn,nbn))
    for l in range(nbe):
        x1=coord[tri[l,0],0]
        x2=coord[tri[l,1],0]
        x3=coord[tri[l,2],0]
        y1=coord[tri[l,0],1]
        y2=coord[tri[l,1],1]
        y3=coord[tri[l,2],1]
        mes1=0.5*abs(((x2-x1)*(y3-y1)-(x3-x1)*(y2-y1)))
        K=coeffelem_p1_rigid(x1,x2,x3,y1,y2,y3,mes1,nbe,coord,tri)
        #la matrice de masse
        ml=masse(x1,x2,x3,y1,y2,y3,mes1,nbe,coord,tri)                          
        #le vecteur de Tl
        fl=coeffelem_p1_source(x1,x2,x3,y1,y2,y3,mes1,nbe,coord,tri)
        I1=tri[l,0]
        I2=tri[l,1]
        I3=tri[l,2]
        A[I1,I1]+=K[0,0]
        A[I2,I2]+=K[1,1]
        A[I3,I3]+=K[2,2]
        A[I2,I1]+=K[1,0]
        A[I3,I1]+=K[2,0]
        A[I1,I2]+=K[0,1]
        A[I3,I2]+=K[2,1]
        A[I1,I3]+=K[0,2]
        A[I2,I3]+=K[1,2]
        F[I1]+=fl[0]
        F[I2]+=fl[1]
        F[I3]+=fl[2] 
        k[I1,I1]=A[I1,I1]
        k[I2,I2]=A[I2,I2]
        k[I3,I3]=A[I3,I3]
        k[I1,I2]=A[I1,I2]
        k[I1,I3]=A[I1,I3]
        k[I2,I1]=A[I2,I1]
        k[I2,I3]=A[I2,I3]
        k[I3,I1]=A[I3,I1]
        k[I3,I2]=A[I3,I2]
    
    for a in range(nba):
        # les sommets du triangle Tl
        x1=coord[ar[a,0],0]
        x2=coord[ar[a,1],0]
        y1=coord[ar[a,0],1]
        y2=coord[ar[a,1],1]
        #la mesure du triangle     
        mes2=sqrt((x1-x2)**2+(y2-y1)**2)
        #le vecteur ea pour l'arret Aa
        ea=coeffelem_P1_transf(x1,x2,y1,y2,mes2,nba,coord,ar)
        # la matrice pa
        Pa=coeffelem_P1_poids(x1,x2,y1,y2,mes2,nba,coord,ar)
        codeFou=refa[a,0]
        if codeFou==3:
            I1=ar[a,0]
            I2=ar[a,1]
            A[I1,I1]+=Pa[0,0]
            A[I2,I2]+=Pa[1,1]
            A[I2,I1]+=Pa[1,0]
            A[I1,I2]+=Pa[0,1]
            F[I1]+=ea[0]
            F[I2]+=ea[1]
#Ajouter les 2 termes de dirichlet"
    for a in range(nba):
        # les sommets du triangle Tl
        x1=coord[ar[a,0],0]
        x2=coord[ar[a,1],0]
        y1=coord[ar[a,0],1]
        y2=coord[ar[a,1],1]
        #la mesure du triangle     
        mes2=sqrt((x1-x2)**2+(y2-y1)**2)
        #le vecteur ea pour l'arret Aa
        ea_dir=coeffelem_P1_transf_dir(x1,x2,y1,y2,mes2,nba,coord,ar)
        # la matrice pa
        Pa_dir=coeffelem_P1_poids_dir(x1,x2,y1,y2,mes2,nba,coord,ar)
        codeDir=refa[a,0] #
        if codeDir==1:
            I1=ar[a,0]
            I2=ar[a,1]
            A[I1,I1]+=Pa_dir[0,0]
            A[I2,I2]+=Pa_dir[1,1]
            A[I2,I1]+=Pa_dir[1,0]
            A[I1,I2]+=Pa_dir[0,1]
            F[I1]+=ea_dir[0]
            F[I2]+=ea_dir[1]
            
    return k,A,F

#resolution: solveur de Cholesky
#Factoristaion
def facto_chol(A):
    n=len(A)
    L=zeros((n,n))
    U=zeros((n,n))
    M=zeros((n,n))
    for i in range(n):
        L[:,i]=A[:,i]/A[i,i]
        U[i,:]=A[i,:]
        for j in range(n):
            for k in range(n):
                M[j,k]=L[j,i]*U[i,k]
        A=A-M
    D=diag(diag(U))
    C=dot(L,sqrt(D))
    return C,C.T

#Descente
def descente(C,F):
    n=len(C)
    y=zeros(n)
    for i in range(n):
        s=0*0
        if i!=0:
            for j in range(i):
                s=s+C[i,j]*y[j]
        y[i]=(F[i]-s)/C[i,i]
    return y

#Remontée
def remontee(U,y):
    n=len(U)
    uh=zeros(n)
    for i in range(n-1,-1,-1):
        s=0*0
        for j in range(n-1,i,-1):
            s=s+U[i,j]*uh[j]
        uh[i]=(y[i]-s)/U[i,i]
    return uh

"""Post-traiteement"""
#Calcul de rh et U
def conv(nbn):
    rh=0.0
    u=zeros(nbn)
    for i in range(nbn):
        psi=zeros(nbn)
        x1=coord[i,0]
        y1=coord[i,1]
        u[i]=fct_u(x1,y1)
        for j in range(nbn):
            #foction chapeau
            if j==i:
                psi[j]=1
                fig=figure(1)
                ax1=fig.add_subplot(111,projection='3d')
                ax1.plot_trisurf(coord[:,0],coord[:,1],tri,psi)
                show()
    return u

##Affichage des resultats numeriques
# affichage du maillage
FichierMaillage = open("./fusee.msh","r")
[nbn,nbe,nba,coord,tri,ar,refn,reft,refa] = charge_et_affiche_maillage(FichierMaillage)
FichierMaillage.close()
pas_et_qualite_maillage(nbe,coord,tri) #(0.013336676905206365, 2.185689303498792)

#affectation de lamatrice EF
[K,A,F]=assemblage_EF_P1(nbn,nbe,nba,tri,ar,refa)
#Factoristaion
[C,Ct]=facto_chol(A)
#Descente
y=descente(C,F)
#Remontée
Uh=remontee(Ct,y) #sol approchée
#convergence
#conv(nbn)

#Graphique
spy(A)
plt.title('')
plt.title('Matrice creuse')
plt.show()

plt.plot(range(nbn),Uh,"r-",label="Uh")
plt.title('')
plt.xlabel('nbn')
plt.ylabel('Uh')
plt.legend()
plt.title('Solution Uh')
plt.show()

