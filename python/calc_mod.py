import scipy.ndimage as nd
import numpy as np
from numpy.linalg import inv
def calc(image,hz,dx):
    if (image[0,0]<128): #Binariza la imagen si es de fondo negro
        imbin=np.where(image<30,0,255)
    else:   #Binariza la imagen si es de fondo blanco
        imbin=np.where(image>200,0,255)
    imfl=nd.median_filter(imbin,(3,3))
    sep,num=nd.label(imfl)
    cm=nd.measurements.center_of_mass(imfl,sep,list(range(1,num+1))) 
    centros=[]    #Lista que contendrá la coordenada "y" de los centros de masa
    for k in range (len(cm)):
        centros.append(cm[k][0])
    centros.pop(0) 
    def min_cuadrados(dt,dy):
        fun=[]
        fun.append(lambda x:np.ones_like(dt))
        fun.append(lambda x:dt)
        fun.append(lambda x:(1/2)*(dt**2))
        Xt=[]
        for fu in fun:
            Xt.append(fu(dt))
        Xt= np.array(Xt)
        X=Xt.transpose()
        a = np.dot(np.dot(inv(np.dot(Xt,X)),Xt),dy)
        return a
    t=np.linspace(1,num-1,num-1)
    y=np.array(centros)
    coef=min_cuadrados(t,y)
    apix=coef[-1]
    acel=apix*(hz**2)*dx #Para tener la aceleración en mm/s^2
    return(acel)
