def calc(image,hz,dx):
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
        return(x)
    def labeli(x): #Cantidad de bolas por arriba
        i=0
        num=[]
        while i<len(x[:,0]):
            if 255 in x[i] and 255 not in x[i-1]:
                num.append(i)
            i+=1
        return(num)
    def labelf(x): #Cantidad de bolas por debajo
        i=0
        num=[]
        while i<len(x[:,0]):
            if 255 not in x[i] and 255 in x[i-1]:
                num.append(i)
            i+=1
        return(num)
    def cm(bi,bf): #Centros de masa
        centros=[[None] for i in range(len(bi)-1)]
        for j in range(len(bi)-1):
            centros[j][0]=((bi[j+1]+bf[j+1])/2)
        return(centros)
    def dot(a,b): #Multiplicación de matrices
        Pr=[[None]*(len(b[0])) for n in range (len(a))]
        cas=0
        for i in range(len(a)):
            for j in range(len(b[1])):
                for k in range(len(a[1])):
                    Pr[i][j]=(a[i][k]*b[k][j])+cas
                    cas=Pr[i][j]
                    # print(Pr[i,j])
                cas=0
        return(Pr)
    def inv(a):   #Inversa de una matriz
        Mt=[[None]*(len(a)) for m in range (len(a[1]))]
        for i in range(len(a[1])):
            for j in range (len(a)):
                Mt[i][j]=a[j][i]            
        det = (a[0][0]*((a[1][1]*a[2][2])-(a[1][2]*a[2][1])))-(a[0][1]*((a[1][0]*a[2][2])-(a[1][2]*a[2][0])))+(a[0][2]*((a[1][0]*a[2][1])-(a[1][1]*a[2][0])))
        Madj=[[None]*len(a) for m in range (len(a[1]))]
        conk=1
        conl=1
        conk2=2
        conl2=2
        for k in range(len(Madj)):
            if k==2:
                conk=-2
                conk2=-1
            if k==1:
                conk2=-1
            for l in range(len(Madj[1])):    
                if l==2:
                    conl=-2
                    conl2=-1
                if l==1:
                    conl2=-1
                Madj[k][l]=((Mt[k+conk][l+conl]*Mt[k+conk2][l+conl2])-(Mt[k+conk][l+conl2]*Mt[k+conl2][l+conl]))/det            
        return(Madj)  
    def min_cuadrados(t,y):
        Xt=[[None]*len(t) for l in range(3)]
        X=[[None]*3 for m in range(len(t))]
        for i in range (len(t)):
            Xt[0][i]=1
            Xt[1][i]=t[i]
            Xt[2][i]=(t[i]**2)/2
            X[i][0]=1
            X[i][1]=t[i]
            X[i][2]=(t[i]**2)/2
        P=dot(Xt,X)
        inversa=inv(P)   
        pr=dot(inversa,Xt)  
        coef=dot(pr,y)
        return(coef)
    binar(image)
    bi=labeli(image)
    bf=labelf(image)
    cen=cm(bi,bf)
    time=[]
    for n in range(len(cen)):
        time.append(n+1)
    apix=min_cuadrados(time,cen)
    acel=apix[2][0]*(hz**2)*dx
    return(acel)
