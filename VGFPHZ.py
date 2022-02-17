from math import pi, cos, sin, e
from numpy import zeros, matmul, loadtxt, array
from matplotlib import pyplot

#needed variables
Ver = []
Theta = []
Hor = []
Phi = []

#value from a sphere on page 42
NDATA = 181 #a parameter that's given

#variables from VGFFMC
V = [0,0,0]
rrr = zeros([3,3],float)
rhat = []#zeros([3],float)
thhat = []#zeros([3],float)
phhat = []#zeros([3],float)
deg = pi/180.0

#values from what Bro. Lines gave me
nuse = 720
w = 10.0 #* 10**-6
alpha = 0.0
beta = 0.0
gamma = 0.0
RAD = 1.55 #*10**-6 #radius given to me
ER = 1.48590004
EI = 0.122000001
temp = 15.6099672
volume = 15.5985308
count = 52

#pull from file
r = loadtxt("figure3fields", skiprows = 5, usecols = (0,1,2), max_rows = 720)
dvol = loadtxt("figure3fields", skiprows = 5, usecols = 3, max_rows = 720)
Ereal = array(loadtxt("figure3fields", skiprows = 725, usecols = (0,2,4), max_rows = 720))
Eimaginary = array(loadtxt("figure3fields", skiprows = 725, usecols = (1,3,5), max_rows = 720))
E = Ereal + 1j * Eimaginary

#rotation matrix
def rot(rrr, alpha, beta, gamma):
    
    rrr[0][0] = cos(alpha)*cos(beta)*cos(gamma)-sin(alpha)*sin(gamma)
    rrr[0][1] = sin(alpha)*cos(beta)*cos(gamma)+cos(alpha)*sin(gamma)
    rrr[0][2] = -sin(beta)*cos(gamma)
    
    rrr[1][0] = -cos(alpha)*cos(beta)*sin(gamma)-sin(alpha)*cos(gamma)
    rrr[1][1] = -sin(alpha)*cos(beta)*sin(gamma)+cos(alpha)*cos(gamma)
    rrr[1][2] = sin(beta)*sin(gamma)
    
    rrr[2][0] = cos(alpha)*sin(beta)
    rrr[2][1] = sin(alpha)*sin(beta)
    rrr[2][2] = cos(beta)
    
    return rrr

#compute chi
eps = complex(ER,EI)
chi = (eps-1)/(4*pi)

#print(dvol)

dsum = 0
for n in range(nuse):
    dsum += dvol[n]
print("Total volume =",dsum)

rot(rrr,alpha,beta,gamma)

#Calculate normalized differential cross section
k = 2 * pi / w

#normalization constant
pa2 = 1/(pi * RAD**2)

#calculate scattering amplitude
PH = 0
#Ingegrate over dOmega
for iang in range(NDATA-1):
    TH = pi*iang/(NDATA-1)

    #make H and V components
    V[0] = cos(TH)*cos(PH)
    V[1] = cos(TH)*sin(PH)
    V[2] = -sin(TH)
    
    thhat = matmul(rrr,V)
    
    V[0] = -sin(PH)
    V[1] = cos(PH)
    V[2] = 0
    
    phhat = matmul(rrr,V)
    
    #make H and V components
    V[0] = sin(TH)*cos(PH)
    V[1] = sin(TH)*sin(PH)
    V[2] = cos(TH)
    
    rhat = matmul(rrr,V)
    
    #print(thhat,phhat)
    
    fh = 0
    fv = 0
    for mu in range(nuse):
        rdot = 0
        
        for j in range(3):
            rdot += r[mu][j]*rhat[j]  
        
        #print(rdot)
        #calculate factors that don't depend on i and j
        c = 1j*k**3*chi*(dvol[mu])*e**(-1j*k*rdot)
        #print(rdot)
        
        #put scattering amplitude together
        tde = 0
        pde = 0
        
        #Calculate thhat dot E and phhat dot E
        for i in range(3):
            tde += thhat[i]*E[mu][i]
            pde += phhat[i]*E[mu][i]
        #print(tde,pde)
        
        #put scattering amplitude together
        fh += c*tde
        fv += c*pde
        #print(fv)
    
    #separate into horizontal and verticle components while converting to dif. scat. cross section
    hor = fh*fh.conjugate() / k**2
    ver = fv*fv.conjugate() / k**2
    #print(hor,ver)
    
    #the complex part will be 0, so we'll make it real now
    hor = hor.real
    ver = ver.real
    
    #normalize and add a point to the plot
    Ver.append(ver*pa2)
    Hor.append(hor*pa2)
    Theta.append(TH/deg)
    Phi.append(PH/deg)

#display results
pyplot.plot(Theta, Ver)
pyplot.title("Vertical vs Theta")
pyplot.xlabel("Scattering angle (degrees)")
pyplot.ylabel("Normalized Scattering Cross Section (x10^-3)")
pyplot.show()
pyplot.plot(Theta, Hor)
pyplot.title("Horizontal vs Theta")
pyplot.xlabel("Scattering angle (degrees)")
pyplot.ylabel("Normalized Scattering Cross Section (x10^-3)")
pyplot.show()

#print(dvol)
