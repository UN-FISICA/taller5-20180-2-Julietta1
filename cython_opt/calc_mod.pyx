from cython.view cimport array as cvarray
from cpython cimport array
import array
def calc(image,hz,dx):
    cdef double [:,:] image_view = image
    def binar(x): #Binariza la matriz
        if x[0,0]<128:
            for i in range(len(x[:,0])):
                for j in range(len(x[0,:])):
                    if x[i,j]<30:
                            x[i,j]=0
                    else:
                            x[i,j]=255
        else:
            for i in range(len(x[:,0])):
                for j in range(len(x[0,:])):
                    if x[i,j]>200:
                            x[i,j]=0
                    else:
                            x[i,j]=255
        cdef double [:,:] x_view = x
        return(x_view)  
    def labelf(x): #Cantidad de bolas por debajo
        cdef int i=0
        cdef int con=0
        cdef array.array a = array.array("i",[])
        cdef array.array num
        num=array.clone(a,len(x[:,0]),0)
        while i<len(x[:,0]):
            if 255 not in x[i] and 255 in x[i-1]:
                num[con]=i
                con+=1
            i+=1
        array.resize(num,con)
        return(num)
    def labeli(x): #Cantidad de bolas por arriba
        cdef int i=0
        cdef int con=0
        cdef array.array a = array.array("i",[])
        cdef array.array num
        num=array.clone(a,len(x[:,0]),0)
        while i<len(x[:,0]):
            if 255 in x[i] and 255 not in x[i-1]:
                num[con]=i-1
                con+=1
            i+=1
        array.resize(num,con)
        return(num)
    def cm(bi,bf): #Centros de masa
        centros=cvarray(shape=(len(bi)-1,1),itemsize=sizeof(float), format="f")
        for j in range(centros.shape[0]):
            centros[j,0]=(bi[j+1]+bf[j+1])/2
        cdef float[:,:] centros_view=centros
        return(centros_view)
    def dot(a,b): #Multiplicación de matrices
        Pr=cvarray(shape=(a.shape[0],b.shape[1]),itemsize=sizeof(float), format="f")
        cdef float cas=0
        for i in range(a.shape[0]):
            for j in range(b.shape[1]):
                for k in range(a.shape[1]):
                    Pr[i,j]=(a[i,k]*b[k,j])+cas
                    cas=Pr[i,j]
                   # print(Pr[i,j])
                cas=0
            #print("fin")
        cdef float [:,:] Pr_view=Pr
        return(Pr_view)
    def inv(a):   #Inversa de una matriz
        Mt=cvarray(shape=(a.shape[1],a.shape[0]),itemsize=sizeof(float), format="f")
        for i in range(a.shape[1]):
            for j in range (a.shape[0]):
                Mt[i,j]=a[j,i]            
        cdef float det = (a[0,0]*((a[1,1]*a[2,2])-(a[1,2]*a[2,1])))-(a[0,1]*((a[1,0]*a[2,2])-(a[1,2]*a[2,0])))+(a[0,2]*((a[1,0]*a[2,1])-(a[1,1]*a[2,0])))
        Madj=cvarray(shape=(a.shape[1],a.shape[0]),itemsize=sizeof(float), format="f")
        conk=1
        conl=1
        conk2=2
        conl2=2
        for k in range(Madj.shape[0]):
            if k==2:
                conk=-2
                conk2=-1
            if k==1:
                conk2=-1
            for l in range(Madj.shape[1]):    
                if l==2:
                    conl=-2
                    conl2=-1
                if l==1:
                    conl2=-1
                Madj[k,l]=((Mt[k+conk,l+conl]*Mt[k+conk2,l+conl2])-(Mt[k+conk,l+conl2]*Mt[k+conl2,l+conl]))/det            
        cdef float [:,:] Mt_view=Mt
        cdef float [:,:] Madj_view=Madj
        return(Madj_view)
    def min_cuadrados(t,y):
        Xt=cvarray(shape=(3,len(t)),itemsize=sizeof(float), format="f")
        X=cvarray(shape=(len(t),3),itemsize=sizeof(float), format="f")    
        for i in range (len(t)):
            Xt[0,i]=1
            Xt[1,i]=t[i]
            Xt[2,i]=(t[i]**2)/2
            X[i,0]=1
            X[i,1]=t[i]
            X[i,2]=(t[i]**2)/2
        P=dot(Xt,X)
        inversa=inv(P)   
        pr=dot(inversa,Xt)  
        coef=dot(pr,y)
        cdef float [:,:] Xt_view= Xt
        cdef float [:,:] X_view=X
        cdef float [:,:] P_view=P
        cdef float [:,:] inversa_view=inversa
        cdef float [:,:] pr_view=pr
        cdef float [:,:] coef_view=coef
        return(coef_view)
    binar(image_view)    
    bf=labelf(image_view)
    bi=labeli(image_view)
    cen=cm(bi,bf)
    cdef array.array t= array.array("f",[])
    cdef array.array time
    time=array.clone(t,cen.shape[0],0)  #Tiempos
    for w in range(len(time)):
        time[w]=w+1
    apix=min_cuadrados(time,cen)
    acel=apix[-1,0]*(hz**2)*dx
    return(acel)
