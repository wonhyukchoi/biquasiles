from sympy import *
import datetime

Z=Symbol('Z')

def mtranspose(M):
    ### transpose matrix### 
    if len(M)==1 and len(M[0])==1: return M
    out=[]
    for i in range(1,len(M[0])+1):
        r = []
        for j in range(1,len(M)+1):
            r.append(M[j-1][i-1])
        out.append(r)
    return out

def avail(M,j,k):
    ### List available entries 
    L = range(1,len(M)+1)
    for i in range(1,len(M)+1):
        if M[j-1][i-1]!=0 and M[j-1][i-1] in L:
            L.remove(M[j-1][i-1])
        if M[i-1][k-1]!=0 and M[i-1][k-1] in L:
            L.remove(M[i-1][k-1])
    return tuple(L)



def tm(M):
    ### tuple matrix
    out  = []
    for x in M: out.append(tuple(x))
    return tuple(out)



def lm(M):
    ### list matrix 
    out = []
    for x in M: out.append(list(x))
    return list(out)


def invn(k,n):
    ### k inverse mod n### 
    out = False
    if n == 0:
        return Rational(1)/k
    for j in range(1,n):
        if j*k %n == 1: out = j
    return out


def getcolumn(M,j):   # get column n from matrix M
    ### Gets column n from matrix M### 
    c =[]
    for i in range(1,len(M)+1):
        c.append(M[i-1][j-1])
    return c


def bfindzero(M):
    ###Find position of first zero entry in a birack matrix###
    for i in range(1,len(M[0])+1):
        for j in range(1,len(M[0])+1):
            if M[0][i-1][j-1] == 0:
                return (0,i,j)
            if M[1][i-1][j-1] == 0:
                 return (1,i,j)
    return False


def bbfindZ(M):
    ### Find position of first zero entry in a birack matrix
    for i in range(1,len(M[0])+1):
        for j in range(1,len(M[0])+1):
            if M[0][i-1][j-1] == Z:
                return (0,i,j)
            if M[1][i-1][j-1] == Z:
                 return (1,i,j)
    return False


def bqslist(X):
    ### Find biquasiles
    working,out=[tuple([tm(X[0]),tm(X[1])])],[]
    while working!=[]:
        w=working[0]
        working[0:1]=[]
        f=bfindzero(w)
        if not f:
            out.append(w)
        if f:
            w2=[lm(w[0]),lm(w[1])]
            for j in avail(w2[f[0]],f[1],f[2]):
                w2[f[0]][f[1]-1][f[2]-1]=j
                w3=bqsfill(w2)
                if w3:
                    working.append(tuple([tm(w3[0]),tm(w3[1])]))
    return out
    

def permtest(L):
    for x in range(0,len(L)):
        for y in range(0,x):
            if L[x]==L[y] and L[x]!=0:
                return False
    return True


def bqsinvcheck(X):
    S,D,St,Dt=X[0],X[1],mtranspose(X[0]),mtranspose(X[1])
    for x in range(0,len(S)):
        if not permtest(S[x]): return False
        if not permtest(St[x]): return False
        if not permtest(D[x]): return False
        if not permtest(Dt[x]): return False
    return True
    


def bqsfill(X):
    ### Fill biquasile matrix with exchange laws
    S,D=lm(X[0]),lm(X[1])
    C=True
    while C:
        C=False
        if not bqsinvcheck(X):
            return False
        for x in range(1,len(S)+1):
            for y in range(1,len(S)+1):
                for a in range(1,len(S)+1):
                    for b in range(1,len(S)+1):
                        if D[a-1][b-1]!=0 and S[y-1][D[a-1][b-1]-1]!=0 and D[x-1][S[y-1][D[a-1][b-1]-1]-1]!=0 and S[a-1][D[x-1][S[y-1][D[a-1][b-1]-1]-1]-1]!=0 and S[a-1][D[x-1][y-1]-1]!=0 and D[S[a-1][D[x-1][y-1]-1]-1][b-1]!=0 and S[y-1][D[S[a-1][D[x-1][y-1]-1]-1][b-1]-1]!=0 and D[x-1][S[y-1][D[S[a-1][D[x-1][y-1]-1]-1][b-1]-1]-1]!=0 and S[S[a-1][D[x-1][y-1]-1]-1][D[x-1][S[y-1][D[S[a-1][D[x-1][y-1]-1]-1][b-1]-1]-1]-1]!=0:
                            if S[a-1][D[x-1][S[y-1][D[a-1][b-1]-1]-1]-1]==0 and S[S[a-1][D[x-1][y-1]-1]-1][D[x-1][S[y-1][D[S[a-1][D[x-1][y-1]-1]-1][b-1]-1]-1]-1]!=0:
                                S[a-1][D[x-1][S[y-1][D[a-1][b-1]-1]-1]-1]=S[S[a-1][D[x-1][y-1]-1]-1][D[x-1][S[y-1][D[S[a-1][D[x-1][y-1]-1]-1][b-1]-1]-1]-1]
                                C=True
                            if S[a-1][D[x-1][S[y-1][D[a-1][b-1]-1]-1]-1]!=0 and S[S[a-1][D[x-1][y-1]-1]-1][D[x-1][S[y-1][D[S[a-1][D[x-1][y-1]-1]-1][b-1]-1]-1]-1]==0:
                                S[S[a-1][D[x-1][y-1]-1]-1][D[x-1][S[y-1][D[S[a-1][D[x-1][y-1]-1]-1][b-1]-1]-1]-1]=S[a-1][D[x-1][S[y-1][D[a-1][b-1]-1]-1]-1]
                                C=True
                            if S[a-1][D[x-1][S[y-1][D[a-1][b-1]-1]-1]-1]!=S[S[a-1][D[x-1][y-1]-1]-1][D[x-1][S[y-1][D[S[a-1][D[x-1][y-1]-1]-1][b-1]-1]-1]-1]:
                                return False
                        if D[x-1][y-1]!=0 and S[a-1][D[x-1][y-1]-1]!=0 and D[S[a-1][D[x-1][y-1]-1]-1][b-1]!=0 and S[y-1][D[S[a-1][D[x-1][y-1]-1]-1][b-1]-1]!=0 and D[a-1][b-1]!=0 and S[y-1][D[a-1][b-1]-1]!=0 and D[x-1][S[y-1][D[a-1][b-1]-1]-1] and S[a-1][D[x-1][S[y-1][D[a-1][b-1]-1]-1]-1] and D[S[a-1][D[x-1][S[y-1][D[a-1][b-1]-1]-1]-1]-1][b-1] and S[S[y-1][D[a-1][b-1]-1]-1][D[S[a-1][D[x-1][S[y-1][D[a-1][b-1]-1]-1]-1]-1][b-1]-1]!=0:
                            if S[y-1][D[S[a-1][D[x-1][y-1]-1]-1][b-1]-1]==0 and S[S[y-1][D[a-1][b-1]-1]-1][D[S[a-1][D[x-1][S[y-1][D[a-1][b-1]-1]-1]-1]-1][b-1]-1]!=0:
                                S[y-1][D[S[a-1][D[x-1][y-1]-1]-1][b-1]-1]=S[S[y-1][D[a-1][b-1]-1]-1][D[S[a-1][D[x-1][S[y-1][D[a-1][b-1]-1]-1]-1]-1][b-1]-1]
                                C=True
                            if S[y-1][D[S[a-1][D[x-1][y-1]-1]-1][b-1]-1]!=0 and S[S[y-1][D[a-1][b-1]-1]-1][D[S[a-1][D[x-1][S[y-1][D[a-1][b-1]-1]-1]-1]-1][b-1]-1]==0:
                                S[S[y-1][D[a-1][b-1]-1]-1][D[S[a-1][D[x-1][S[y-1][D[a-1][b-1]-1]-1]-1]-1][b-1]-1]=S[y-1][D[S[a-1][D[x-1][y-1]-1]-1][b-1]-1]
                                C=True
                            if S[y-1][D[S[a-1][D[x-1][y-1]-1]-1][b-1]-1]!=S[S[y-1][D[a-1][b-1]-1]-1][D[S[a-1][D[x-1][S[y-1][D[a-1][b-1]-1]-1]-1]-1][b-1]-1]:
                                return False
    return tuple([tm(S),tm(D)])


def reptest(p):   
    ### Test whether p has repeated non-zero entries ###
    q = True
    L = []
    for i in range(1,len(p)+1):
        if p[i-1] != 0:
            if p[i-1] in L:
                q = False
            else:
                L.append(p[i-1])
    return q



def bisofill(M,N,f):
    ### fill isomorphism f:M -> N ###
    m,n = len(M[0]), len(N[0])
    if m != n: return False
    keepgoing = True
    while keepgoing:
        keepgoing = False
        for i in range(1,n+1):
            for j in range(1,n+1):
                if f[i-1] !=0 and f[j-1] != 0:
                    if f[M[0][i-1][j-1]-1] == 0: 
                        f[M[0][i-1][j-1]-1] = N[0][f[i-1]-1][f[j-1]-1]
                        keepgoing = True
                    if N[0][f[i-1]-1][f[j-1]-1] != f[M[0][i-1][j-1]-1]: return False
                    if f[M[1][i-1][j-1]-1] == 0: 
                        f[M[1][i-1][j-1]-1] = N[1][f[i-1]-1][f[j-1]-1]
                        keepgoing = True
                    if N[1][f[i-1]-1][f[j-1]-1] != f[M[1][i-1][j-1]-1]: return False
    if reptest(f): return f
    return False


def bhomfill(M,N,f):
    ### fill homomorphism f:M -> N ###
    m,n = len(M[0]), len(N[0])
#    if m != n: return False
    keepgoing = True
    while keepgoing:
        keepgoing = False
        for i in range(1,m+1):
            for j in range(1,m+1):
                if f[i-1] !=0 and f[j-1] != 0:
                    if f[M[0][i-1][j-1]-1] == 0: 
                        f[M[0][i-1][j-1]-1] = N[0][f[i-1]-1][f[j-1]-1]
                        keepgoing = True
                    if N[0][f[i-1]-1][f[j-1]-1] != f[M[0][i-1][j-1]-1]: return False
                    if f[M[1][i-1][j-1]-1] == 0: 
                        f[M[1][i-1][j-1]-1] = N[1][f[i-1]-1][f[j-1]-1]
                        keepgoing = True
                    if N[1][f[i-1]-1][f[j-1]-1] != f[M[1][i-1][j-1]-1]: return False
    return f



def bhomlist(M,N):
    ### find biquandle homomorphisms ###
    z = []
    for i in range(1,len(M[0])+1): z.append(0)
    L,out = [z],[]
    while len(L) > 0:
        w = L[0]
        L[0:1] = []  
        i = hfindzero(w)
        if not i: out.append(list(w))
        else:    
            for j in range(1,len(N[0])+1):
                phi = list(w)
                phi[i-1] = j
                v = bhomfill(M,N,phi)
                if v: L.append(tuple(v))
    return out

     

def bisotest(M,N):
    ###test for biquandle isomorphism###
    z = []
    for i in range(1,len(M[0])+1): z.append(0)
    L,out = [z],[]
    while len(L) != 0:
        w = L[0]
        L[0:1] = []  
        i = hfindzero(w)
        if (not i) and permtest(w): return True
        else:
            for j in pavail(w):
                phi = list(w)
                phi[i-1] = j
                v = bisofill(M,N,phi)
                if v: L.append(tuple(v))
    return False


def breducelist(L):
    ###Remove isomorphic copies from L###
    out = [tm(L[0])]
    W = L
    while len(W)>0:
        x = W[0]
        W[0:1] = []
        newbiq = True
        for y in out:
            if bisotest(x,y): newbiq = False
        if newbiq: out.append(tm(x))
    return out

def hfindzero(f):   
    ### find zero in homomorphism template ###
    j = -1
    for i in range(0,len(f)):
        if f[i] == 0: 
            j = i+1
            break
    if j < 0: out = False
    else: out = j
    return out

def pavail(v):
    ### List available entries ###
    L = []
    for i in range(1,len(v)+1):
        if not i in v: L.append(i)
    return tuple(L)


def pdhfindzero(x):
    for j in range(0,2):
        for k in range(0,len(x[0])):
            if x[j][k]==0:
                return([j,k])
    return False

def pdcolors(X,D):
    ### find biquasile colorings of D
    working,out,x=[],[],[[],[]]
    for j in range(0,2*len(D)):
        x[0].append(0)
        x[1].append(0)
    working=[tm(x)]
    while working!=[]:
        w=working[0]
        working[0:1]=[]
        f=pdhfindzero(w)
        if not f:
            if pdcolorfill(X,D,w):
                out.append(w)
        else:
            for k in range(1,len(X[0])+1):
                w2=lm(w)
                w2[f[0]][f[1]]=k
                w3=pdcolorfill(X,D,w2)
                if w3:
                    working.append(tm(w3))
    return out
        
def bqknotlist(X):
	### Finds biquasile colorings for each knot
	### args: 
	### X = biquasile 
	### returns:
	### L = matrix with knots and their coloring numbers 
	L = []
	for x in gknotlist: 
		L.append([x, len(pdcolors(X,gauss2pd(gknot[x])))])
	return L

def bqslinklist(X):
	### Finds biquasile colorings for each link
	### args: 
	### X = biquasile 
	### returns:
	### L = matrix with links and their coloring numbers 
	L = []
	for x in glinklist: 
		L.append([x, len(pdcolors(X,gauss2pd(glink[x])))])
	return L
def pdcolorfill(X,G,w1):
    ### fill coloring 
    S,D=X[0],X[1]
    C=True
    w=lm(w1)
    while C:
        C=False
        for x in G:
            if x[0]==1:
                if w[0][x[1]-1] !=0 and w[0][x[2]-1]==0:
                    w[0][x[2]-1]=w[0][x[1]-1]
                    C=True
                if w[0][x[1]-1] ==0 and w[0][x[2]-1]!=0:
                    w[0][x[1]-1]=w[0][x[2]-1]
                    C=True
                if w[0][x[1]-1]!=w[0][x[2]-1]:
                    return False
                if w[0][x[4]-1] ==0 and w[1][x[1]-1]!=0:
                    w[0][x[4]-1]=w[1][x[1]-1]
                    C=True
                if w[0][x[4]-1] !=0 and w[1][x[1]-1]==0:
                    w[1][x[1]-1]=w[0][x[4]-1]
                    C=True
                if w[1][x[1]-1] != w[0][x[4]-1]:
                    return False
                if w[1][x[2]-1] !=0 and w[0][x[3]-1]==0:
                    w[0][x[3]-1]=w[1][x[2]-1]
                    C=True
                if w[1][x[2]-1] ==0 and w[0][x[3]-1]!=0:
                    w[1][x[2]-1]=w[0][x[3]-1]
                    C=True
                if w[1][x[2]-1]!=w[0][x[3]-1]:
                    return False
                if w[1][x[4]-1] !=0 and w[1][x[3]-1]==0:
                    w[1][x[3]-1]=w[1][x[4]-1]
                    C=True
                if w[1][x[4]-1] ==0 and w[1][x[3]-1]!=0:
                    w[1][x[4]-1]=w[1][x[3]-1]
                    C=True
                if w[1][x[4]-1]!=w[1][x[3]-1]:
                    return False
                if w[0][x[1]-1]!=0 and w[0][x[4]-1]!=0 and w[1][x[4]-1]!=0:
                    if w[0][x[3]-1]==0:
                        w[0][x[3]-1]=S[w[0][x[4]-1]-1][D[w[0][x[1]-1]-1][w[1][x[4]-1]-1]-1]
                        C=True
                    if w[0][x[3]-1]!=S[w[0][x[4]-1]-1][D[w[0][x[1]-1]-1][w[1][x[4]-1]-1]-1]:
                        return False
            if x[0]==-1:
                if w[0][x[2]-1]==0 and w[0][x[3]-1]!=0:
                    w[0][x[2]-1]=w[0][x[3]-1]
                    C=True
                if w[0][x[2]-1]!=0 and w[0][x[3]-1]==0:
                    w[0][x[3]-1]=w[0][x[2]-1]
                    C=True
                if w[0][x[2]-1]!=w[0][x[3]-1]:
                    return False
                if w[1][x[2]-1]==0 and w[0][x[1]-1]!=0:
                    w[1][x[2]-1]=w[0][x[1]-1]
                    C=True
                if w[1][x[2]-1]!=0 and w[0][x[1]-1]==0:
                    w[0][x[1]-1]=w[1][x[2]-1]
                    C=True
                if w[1][x[2]-1]!=w[0][x[1]-1]:
                    return False
                if w[1][x[3]-1]==0 and w[0][x[4]-1]!=0:
                    w[1][x[3]-1]=w[0][x[4]-1]
                    C=True
                if w[1][x[3]-1]!=0 and w[0][x[4]-1]==0:
                    w[0][x[4]-1]=w[1][x[3]-1]
                    C=True
                if w[1][x[3]-1]!=w[0][x[4]-1]:
                    return False
                if w[1][x[1]-1]==0 and w[1][x[4]-1]!=0:
                    w[1][x[1]-1]=w[1][x[4]-1]
                    C=True
                if w[1][x[1]-1]!=0 and w[1][x[4]-1]==0:
                    w[1][x[4]-1]=w[1][x[1]-1]
                    C=True
                if w[1][x[1]-1]!=w[1][x[4]-1]:
                    return False
                if w[0][x[3]-1]!=0 and w[0][x[4]-1]!=0 and w[1][x[4]-1]!=0:
                    if w[0][x[1]-1]==0:
                        w[0][x[1]-1]=S[w[0][x[4]-1]-1][D[w[0][x[3]-1]-1][w[1][x[4]-1]-1]-1]
                        C=True
                    if w[0][x[1]-1]!=S[w[0][x[4]-1]-1][D[w[0][x[3]-1]-1][w[1][x[4]-1]-1]-1]:
                        return False
    return tm(w)


gknot = {
  (0):[[-1,1]],
(3,1):[[-1,2,-3,1,-2,3]],
(4,1):[[-1.5,2.5,-3,4,-2.5,1.5,-4,3]],
(5,1):[[-1,2,-3,4,-5,1,-2,3,-4,5]],
(5,2):[[-1,2,-3,4,-5,3,-2,1,-4,5]],
(6,1):[[-1,2,-3.5,4.5,-5.5,6.5,-2,1,-6.5,5.5,-4.5,3.5]],
(6,2):[[-1,2.5,-3.5,4.5,-5.5,1,-6,3.5,-4.5,5.5,-2.5,6]],
(6,3):[[1,-2.5,3.5,-4.5,2.5,-5,6,-1,5,-3.5,4.5,-6]],
#(6,3):[[-1,2,-3.5,4.5,-5.5,1,-6.5,3.5,-4.5,6.5,-2,5.5]],
(6,4):[[-1,2,-3,1,-2,3,-4,5,-6,4,-5,6]],
(6,5):[[-1,2,-3,1,-2,3,4.5,-5.5,6.5,-4.5,5.5,-6.5]],
(7,1):[[-1,2,-3,4,-5,6,-7,1,-2,3,-4,5,-6,7]],
(7,2):[[-1.5,6.5,-7.5,1.5,-2.5,3.5,-4.5,5.5,-6.5,7.5,-5.5,4.5,-3.5,2.5]],
(7,3):[[-1,2,-3,4,-7,6,-5,1,-2,3,-4,5,-6,7]],
(7,4):[[-1,2,-3,4,-5,6,-7,1,-4,3,-2,7,-6,5]],
(7,5):[[-1.5,2.5,-3.5,4.5,-5.5,1.5,-2.5,3.5,-6.5,7.5,-4.5,5.5,-7.5,6.5]],
(7,6):[[-1.5,2.5,-3,4,-5.5,6.5,-4,3,-7.5,1.5,-6.5,5.5,-2.5,7.5]],
(7,7):[[-1.5,2.5,-3,4,-2.5,5.5,-6,7,-5.5,1.5,-4,3,-7,6]],
(8,1):[[-1.5,2.5,-3,4,-2.5,1.5,-5.5,6.5,-7.5,8.5,-4,3,-8.5,7.5,-6.5,5.5]],
(8,2):[[1,-2,3.5,-4.5,5.5,-6.5,7.5,-8.5,2,-1,8.5,-3.5,4.5,-5.5,6.5,-7.5]],
(8,3):[[-1,2,-3,4,-5.5,6.5,-7.5,8.5,-4,3,-2,1,-8.5,7.5,-6.5,5.5]],
(8,4):[[-1.5,2,-3,4,-5,1.5,-6.5,7.5,-8.5,5,-4,3,-2,6.5,-7.5,8.5]],
(8,5):[[-1.5,2.5,-3,4,-5,6,-7,8,-2.5,1.5,-6,7,-8,3,-4,5]],
(8,6):[[1,-2,3.5,-4.5,5.5,-6.5,7.5,-8.5,2,-1,8.5,-7.5,6.5,-3.5,4.5,-5.5]],
(8,7):[[-1.5,2.5,-3.5,1.5,-4,5,-6,7,-8,4,-2.5,3.5,-5,6,-7,8]],
(8,8):[[1,-2,3,-4.5,5.5,-6.5,4.5,-7,8,-5.5,6.5,-3,2,-1,7,-8]],
(8,9):[[-1.5,2.5,-3.5,4.5,-5,6,-7,8,-2.5,3.5,-4.5,1.5,-8,5,-6,7]],
#(8,10):[[-1,2,-3,4,-5,6,-7,1,-2,8.5,-6,7,-8.5,3,-4,5]],
(8,10):[[-1,2,-3,4,-5,1,-2,6.5,-7.5,3,-4,8.5,-6.5,7.5,-8.5,5]],
(8,11):[[1,-2,3.5,-4.5,5.5,-6.5,7.5,-8.5,2,-1,8.5,-5.5,4.5,-3.5,6.5,-7.5]],
(8,12):[[-1,2,-3,4,-5.5,6.5,-4,3,-7.5,8.5,-2,1,-8.5,7.5,-6.5,5.5]],
(8,13):[[1.5,-2.5,3,-4,5,-6,7,-3,8.5,-1.5,2.5,-8.5,4,-7,6,-5]],
(8,14):[[-1,2,-3.5,4.5,-5.5,6.5,-7.5,5.5,-4.5,8.5,-2,1,-6.5,7.5,-8.5,3.5]],
(8,15):[[1.5,-2.5,3.5,-4.5,5.5,-3.5,6.5,-7.5,8.5,-6.5,2.5,-1.5,7.5,-8.5,4.5,-5.5]],
(8,16):[[1.5,-2,3,-4.5,5.5,-6,2,-7.5,4.5,-5.5,8.5,-1.5,7.5,-3,6,-8.5]],
###(8,17):[[1.5,-2.5,3,-4,5.5,-1.5,6.5,-3,4,-7.5,8.5,-5.5,2.5,-6.5,7.5,-8.5]],
(8,17):[[-1.5,2.5,-3,4,-5,6,-2.5,7.5,-4,5,-8.5,1.5,-6,3,-7.5,8.5]],
(8,18):[[-1,2.5,-3.5,4,-5,6.5,-2.5,7,-4,8.5,-6.5,1,-7,3.5,-8.5,5]],
(8,19):[[-1,2,-3,-4,5,1,-2,-6,4,7,-8,-5,6,3,-7,8]],
(8,20):[[-1,2.5,3,-4,-5.5,1,6.5,-3,4,-7.5,8.5,5.5,-2.5,-6.5,7.5,-8.5]],
(8,21):[[-1,2.5,-3.5,-4.5,5.5,-6,-7.5,1,4.5,-5.5,8.5,7.5,-2.5,3.5,6,-8.5]],
(9,2):[[-1.5,2.5,-3.5,4.5,-5.5,6.5,-7.5,8.5,-9.5,1.5,-2.5,9.5,-8.5,7.5,-6.5,5.5,-4.5,3.5]],
(9,24):[[-1.5,2.5,-3,4.5,-5.5,6,-7,8,-2.5,1.5,-6,7,-8,9.5,-4.5,5.5,-9.5,3]],
(9,32):[[-1.5,2.5,-3,4,-5,6,-7,8,-2.5,9.5,-4,7,-8,3,-9.5,1.5,-6,5]],
(10,132):[[1.5,-2.5,3.5,-1.5,-4,5.5,6,-7.5,-8,4,9.5,-6,-10.5,8,2.5,-3.5,7.5,10.5,-5.5,-9.5]],
###(11,34):[[-1,2,-3,4.5,-5.5,6.5,-7.5,1,-8.5,9.5,-2,3,-10.5,8.5,-9.5,10.5,-11,7.5,-4.5,5.5,-6.5,11]],
(11,1):[[-1,2.5,-3.5,-4,5,6.5,-2.5,7,-8,3.5,-6.5,1,-7,8,9.5,-10.5,11.5,-5,4,-9.5,10.5,-11.5]],
(11,2):[[-1,2.5,-3.5,-4.5,5.5,-6.5,7,-8,4.5,-5.5,6.5,9.5,-2.5,10,-11,3.5,-9.5,1,-10,11,8,-7]],
###(11,34):[[-1,2,-3,4.5,5.5,-6.5,7.5,-5.5,8,-9,6.5,-7.5,-10.5,3,-11,-8,9,1,-2,10.5,-4.5,11]],
(11,42):[[-1,2,3.5,-4,5,-6,7,-3.5,8.5,-5,6,9.5,-10.5,11.5,-9.5,1,-2,10.5,-11.5,-7,4,-8.5]]
}

gknotlist=((0),(3,1),(4,1),(5,1),(5,2),(6,1),(6,2),(6,3),(7,1),(7,2),(7,3),(7,4),(7,5),(7,6),(7,7),
	(8,1),(8,2),(8,3),(8,4),(8,5),(8,6),(8,7),(8,8),(8,9),(8,10),(8,11),(8,12),(8,13),(8,14),(8,15),
	(8,16),(8,17),(8,18),(8,19),(8,20),(8,21),(9,32), (10,132), (11,1), (11,2),(11,42))

glinklist=((2,0,1),(4,0,1),(5,0,1),(6,0,1),(6,0,2),(6,0,3),(6,0,4),(6,0,5),(6,1,1),(7,0,1),(7,0,2),(7,0,3),(7,0,4),(7,0,5),(7,0,6),(7,0,7),(7,1,1),(7,1,2))

glink = {
(0,2):[[-1,1],[-2,2]],
(0,3):[[-1,1],[-2,2],[-3,3]],
(2,0,1):[[-1,2],[-2,1]],
(4,0,1):[[-1,2,-3,4],[-4,3,-2,1]],
(5,0,1):[[-1,2.5,-3.5,4],[1,-4,5.5,-2.5,3.5,-5.5]],
(6,0,1):[[-1,2,-3,4],[1,-5.5,6.5,-4,3,-6.5,5.5,-2]],
#(6,0,2):[[-1,2,-5.5,4,-3.5,6],[1,-2,3.5,-4,5.5,-6]],
(6,0,2):[[-1,2,-3,4,-5,6],[1,-2,3,-6,5,-4]],
(6,0,3):[[-1,2,-3,4,-5,6],[1,-2,3,-4,5,-6]],
(6,0,4):[[-1,2,-3,4],[-5,1,-6,3],[-2,6,-4,5]],
#(6,0,5):[[-1,2.5,-3.5,4],[1,-5,6,-2.5],[3.5,-6,5,-4]],
(6,0,5):[[-1,2,-3,4],[1,-5,6,-2],[3,-6,5,-4]],
(6,1,1):[[-1,-2,3,4],[1,-5,-3,6],[2,5,-4,-6]],
(7,0,1):[[-1.5,2.5,-3,4,-5.5,1.5,-6,3,-7.5,5.5],[-2.5,7.5,-4,6]],
(7,0,2):[[-1,2,-3,4,-5,1,-2,5,-6,7],[3,-7,6,-4]],
(7,0,3):[[-1.5,2.5,-3.5,4.5,-5.5,1.5,-2.5,3.5,-6,7],[-4.5,5.5,-7,6]],
(7,0,4):[[-1.5,2.5,-3.5,4.5,-5.5,3.5,-2.5,1.5,-6,7],[-4.5,5.5,-7,6]],
(7,0,5):[[-1,2.5,-3.5,1,-4,5,-6,7],[-2.5,4,-7,6,-5,3.5]],
(7,0,6):[[-1,2.5,-3.5,4.5,-2.5,5,-6,7],[1,-7,6,-5,3.5,-4.5]],
(7,0,7):[[-1,2,-3,4],[1,-5.5,6.5,-2],[3,-7.5,5.5,-6.5,7.5,-4]],
(7,1,1):[[-1,2.5,3,-4.5,-5,1,6.5,-3,-7.5,5],[-2.5,-6.5,4.5,7.5]],
(7,1,2):[[-1,2.5,3,-4.5,-5,1,-6,-3,7,5],[-2.5,6,4.5,-7]]
#2.4:[[-1,4,-5,6,-3,2],[1,-2,3,-4,5,-6]],
#2.5:[[-1,2,-3,4,-5,6,-7,8],[1,-2,3,-8,7,-4,5,-6]],
#3.1:[[-1,2],[-2,1],[-3,3]],
#3.2:[[-1,2],[-2,3,-4,1],[-3,4]],
}

def refgauss(G):
    ### Gauss code of mirror image ###
    out=[]
    for x in G:
        comp=[]
        for y in x: 
            if (y % 1) == 0 and (y>0):
                comp.append(y+0.5)
            if (y % 1) == 0 and (y<0):
                comp.append(y-0.5)
            if (y % 1) == 0.5 and (y>0):
                comp.append(y-0.5)
            if (y % 1) == 0.5 and (y<0):
                comp.append(y+0.5)
        out.append(comp)
    return out
 
 
def revgauss(G):
    ### Gauss code of reverse ###
    out=[]
    for x in G:
        comp=[]
        for i in range(0,len(x)):
            comp.append(x[len(x)-i-1])
        out.append(comp)
    return out     

def transformGC(G):
	###Transforms knot into an array of its mirror image and reverse
	### args: 
	### G = the knot
	### returns:
	### An array of knots
	G = gknot[G]
	return [G, refgauss(G), revgauss(G),revgauss(refgauss(G))]

def transformGC2(G):
	###Transforms links into an array of its mirror image and reverse
	### args: 
	### G = the link
	### returns:
	### An array of link
	G = glink[G]
	return [G, refgauss(G), revgauss(G),revgauss(refgauss(G))]	

def gauss2pd(G):
    ###convert Gauss code to planar diagram###
    out = []
    semiarccount = 0
    complen = []
    for i in range(0,len(G)):
        semiarccount = semiarccount + len(G[i])
        complen.append(len(G[i]))
    nx, cr = [], [] 
    for k in range(0,semiarccount):
        nx.append(0)
        cr.append(0)
    current = 1
    i = 0
    for x in range(1,len(G)+1):
        for y in range(1,len(G[x-1])+1):
            i = i+1
            if y == complen[x-1]:
                nx[i-1] = current
            else:
                nx[i-1] = i+1
            j = 0
            for z in range(1,len(G)+1):
                for w in range(1,len(G[z-1])+1):
                    j = j +1
                    if G[z-1][w-1] + G[x-1][y-1] == 0:
                        cr[i-1] = j 
        current = current + complen[x-1]
    i = 0
    for x in range(1,len(G)+1):
        for y in range(1,len(G[x-1])+1):
            i = i+1
            if G[x-1][y-1] < 0:
                if G[x-1][y-1] % 1 == 0:
                    out.append([ 1,i,nx[cr[i-1]-1],nx[i-1],cr[i-1] ])
                else:
                    out.append([-1,i,cr[i-1],nx[i-1],nx[cr[i-1]-1] ])
    return out



def pdmirror(PD):
    ### Mirror PD
    out=[]
    for x in PD:
        if x[0]==1:
            out.append([-1,x[4],x[1],x[2],x[3]])
        if x[0]==-1:
            out.append([1,x[2],x[3],x[4],x[1]])
    return out

def rdiv(X,x,y):
    ### return x/y
    for z in range(1,len(X)+1):
        if X[z-1][y-1]==x:
            return z
    return False

def ldiv(X,x,y):
    ### return x\y
    for z in range(1,len(X)+1):
        if X[x-1][z-1]==y:
            return z
    return False


def rdmatrix(X):
    ### right division matrix
    out=[]
    for x in range(1,len(X)+1):
        row=[]
        for y in range(1,len(X)+1):
            row.append(rdiv(X,x,y))
        out.append(row)
    return out


def ldmatrix(X):
    ### left division matrix
    out=[]
    for x in range(1,len(X)+1):
        row=[]
        for y in range(1,len(X)+1):
            row.append(ldiv(X,y,x))
        out.append(row)
    return out

def lm2(M):
    ### list matrix### 
    out = []
    for x in M: out.append(lm(x))
    return list(out)


def tm2(M):
    ### list matrix### 
    out = []
    for x in M: out.append(tm(x))
    return tuple(out)



def zmatrix(n):
    ### return n by n Z matrix### 
    out = []
    for i in range(0,n): 
        r = []
        for j in range(0,n): r.append(Z)
        out.append(r)
    return out

def findZ(M):
    ### Find position of first Z entry in a cocycle matrix
    for i in range(1,len(M[0])+1):
        for j in range(1,len(M[0])+1):
            for k in range(1,len(M[0])+1):
                if M[i-1][j-1][k-1] == Z:
                    return (i,j,k)
    return False



def bwlist(X,n):
    ### find Boltzmann weights over Z_n
    S,D,working,out=X[0],X[1],[],[]
    SRD,SLD,DRD,DLD=rdmatrix(S),ldmatrix(S),rdmatrix(D),ldmatrix(D)
    temp=[]
    for j in range(1,len(S)+1):
        temp.append(lm(zmatrix(len(S))))
    for x in range(1,len(S)+1):
        for a in range(1,len(S)+1):
            temp[x-1][a-1][ldiv(D,a,ldiv(S,x,x))-1]=0
            temp[x-1][rdiv(D,ldiv(S,x,x),a)-1][a-1]=0
    working=[tm2(temp)]
    while working!=[]:
        w=working[0]
        working[0:1]=[]
        f=findZ(w)
        if not f:
            if bwfill(X,n,w):
                out.append(tm2(w))
        else: 
            for k2 in range(0,n):
                w2=lm2(w)
                w2[f[0]-1][f[1]-1][f[2]-1]=k2
                w3=bwfill(X,n,w2)
                if w3:
                    working.append(tm2(w3))
    return out



def bwlistN(X,n,N):
    ### find N Boltzmann weights over Z_n
    S,D,working,out=X[0],X[1],[],[]
    SRD,SLD,DRD,DLD=rdmatrix(S),ldmatrix(S),rdmatrix(D),ldmatrix(D)
    temp=[]
    for j in range(1,len(S)+1):
        temp.append(lm(zmatrix(len(S))))
    for x in range(1,len(S)+1):
        for a in range(1,len(S)+1):
            temp[x-1][a-1][ldiv(D,a,ldiv(S,x,x))-1]=0
            temp[x-1][rdiv(D,ldiv(S,x,x),a)-1][a-1]=0
    working=[tm2(temp)]
    while working!=[] and len(out)<N:
        w=working[0]
        working[0:1]=[]
        f=findZ(w)
        if not f:
            if bwfill(X,n,w):
                out.append(tm2(w))
        else: 
            for k2 in range(0,n):
                w2=lm2(w)
                w2[f[0]-1][f[1]-1][f[2]-1]=k2
                w3=bwfill(X,n,w2)
                if w3:
                    working.append(tm2(w3))
    return out

def bwlistTest(X,n):
	### modifying... 
    S,D,working,out=X[0],X[1],[],[]
    SRD,SLD,DRD,DLD=rdmatrix(S),ldmatrix(S),rdmatrix(D),ldmatrix(D)
    temp=[]
    for j in range(1,len(S)+1):
        temp.append(lm(zmatrix(len(S))))
    for x in range(1,len(S)+1):
        for a in range(1,len(S)+1):
            temp[x-1][a-1][ldiv(D,a,ldiv(S,x,x))-1]=0
            temp[x-1][rdiv(D,ldiv(S,x,x),a)-1][a-1]=0
    working=[tm2(temp)]
    return working 
    before = datetime.datetime.now()
    print("We started at:" + str(before))
    o = 0
    while working!=[]:
        w=working[0]
        working[0:1]=[]
        f=findZ(w)      
        if not f:
        	if o == 0:
        		middle = datetime.datetime.now()
        		print("We append at:" + str(middle))
        		o = o +1
        	if bwfill(X,n,w):
        		out.append(tm2(w))	
        else: 
            for k2 in range(0,n):
                w2=lm2(w)
                w2[f[0]-1][f[1]-1][f[2]-1]=k2
                w3=bwfill(X,n,w2)
                if w3:
                    working.append(tm2(w3))
    after = datetime.datetime.now()                
    print("We stop at:" + str(after))
    time1 = middle.day*1440 + middle.hour*60 + middle.minute - before.day*1440 - before.hour*60 - before.minute 
    time2 = after.day*1440 + after.hour*60 + after.minute - middle.day * 1440 - middle.hour*60 - middle.minute
    print("Total time:" + str(time1 + time2) + "\n" +
    	"first half:" + str(time1) + "second half:" +str(time2))                
    return out

def bwfill(X,n,T):
    ### fill with Boltzmann weight conditions
    S,D,M,SRD,SLD,DRD,DLD=X[0],X[1],lm2(T),rdmatrix(X[0]),ldmatrix(X[0]),rdmatrix(X[1]),ldmatrix(X[1])
    keepgoing=True
    while keepgoing:
        keepgoing=False
        for x in range(1,len(S)+1):
            for a in range(1,len(S)+1):
                for y in range(1,len(S)+1):    
                    for b in range(1,len(S)+1):    
                        if M[x-1][a-1][b-1]==Z and M[b-1][S[x-1][D[a-1][b-1]-1]-1][y-1]!=Z and M[S[x-1][D[a-1][b-1]-1]-1][a-1][S[b-1][D[S[x-1][D[a-1][b-1]-1]-1][y-1]-1]-1]!=Z and M[b-1][x-1][y-1]!=Z and M[x-1][a-1][S[b-1][D[x-1][y-1]-1]-1]!=Z and M[S[b-1][D[x-1][y-1]-1]-1][S[x-1][D[a-1][S[b-1][D[x-1][y-1]-1]-1]-1]-1][y-1]!=Z:
                            M[x-1][a-1][b-1]= (-M[b-1][S[x-1][D[a-1][b-1]-1]-1][y-1] -M[S[x-1][D[a-1][b-1]-1]-1][a-1][S[b-1][D[S[x-1][D[a-1][b-1]-1]-1][y-1]-1]-1] +M[b-1][x-1][y-1] +M[x-1][a-1][S[b-1][D[x-1][y-1]-1]-1]+M[S[b-1][D[x-1][y-1]-1]-1][S[x-1][D[a-1][S[b-1][D[x-1][y-1]-1]-1]-1]-1][y-1])%n
                            keepgoing=True
                        if M[x-1][a-1][b-1]!=Z and M[b-1][S[x-1][D[a-1][b-1]-1]-1][y-1]==Z and M[S[x-1][D[a-1][b-1]-1]-1][a-1][S[b-1][D[S[x-1][D[a-1][b-1]-1]-1][y-1]-1]-1]!=Z and M[b-1][x-1][y-1]!=Z and M[x-1][a-1][S[b-1][D[x-1][y-1]-1]-1]!=Z and M[S[b-1][D[x-1][y-1]-1]-1][S[x-1][D[a-1][S[b-1][D[x-1][y-1]-1]-1]-1]-1][y-1]!=Z:
                            M[b-1][S[x-1][D[a-1][b-1]-1]-1][y-1]= (-M[x-1][a-1][b-1] -M[S[x-1][D[a-1][b-1]-1]-1][a-1][S[b-1][D[S[x-1][D[a-1][b-1]-1]-1][y-1]-1]-1] +M[b-1][x-1][y-1] +M[x-1][a-1][S[b-1][D[x-1][y-1]-1]-1]+M[S[b-1][D[x-1][y-1]-1]-1][S[x-1][D[a-1][S[b-1][D[x-1][y-1]-1]-1]-1]-1][y-1])%n
                            keepgoing=True
                        if M[x-1][a-1][b-1]!=Z and M[b-1][S[x-1][D[a-1][b-1]-1]-1][y-1]!=Z and M[S[x-1][D[a-1][b-1]-1]-1][a-1][S[b-1][D[S[x-1][D[a-1][b-1]-1]-1][y-1]-1]-1]==Z and M[b-1][x-1][y-1]!=Z and M[x-1][a-1][S[b-1][D[x-1][y-1]-1]-1]!=Z and M[S[b-1][D[x-1][y-1]-1]-1][S[x-1][D[a-1][S[b-1][D[x-1][y-1]-1]-1]-1]-1][y-1]!=Z:
                            M[S[x-1][D[a-1][b-1]-1]-1][a-1][S[b-1][D[S[x-1][D[a-1][b-1]-1]-1][y-1]-1]-1]= (-M[b-1][S[x-1][D[a-1][b-1]-1]-1][y-1] -M[x-1][a-1][b-1] +M[b-1][x-1][y-1] +M[x-1][a-1][S[b-1][D[x-1][y-1]-1]-1]+M[S[b-1][D[x-1][y-1]-1]-1][S[x-1][D[a-1][S[b-1][D[x-1][y-1]-1]-1]-1]-1][y-1])%n
                            keepgoing=True
                        if M[x-1][a-1][b-1]!=Z and M[b-1][S[x-1][D[a-1][b-1]-1]-1][y-1]!=Z and M[S[x-1][D[a-1][b-1]-1]-1][a-1][S[b-1][D[S[x-1][D[a-1][b-1]-1]-1][y-1]-1]-1]!=Z and M[b-1][x-1][y-1]==Z and M[x-1][a-1][S[b-1][D[x-1][y-1]-1]-1]!=Z and M[S[b-1][D[x-1][y-1]-1]-1][S[x-1][D[a-1][S[b-1][D[x-1][y-1]-1]-1]-1]-1][y-1]!=Z:
                            M[b-1][x-1][y-1]= (M[b-1][S[x-1][D[a-1][b-1]-1]-1][y-1] +M[x-1][a-1][b-1] +M[S[x-1][D[a-1][b-1]-1]-1][a-1][S[b-1][D[S[x-1][D[a-1][b-1]-1]-1][y-1]-1]-1] -M[x-1][a-1][S[b-1][D[x-1][y-1]-1]-1]-M[S[b-1][D[x-1][y-1]-1]-1][S[x-1][D[a-1][S[b-1][D[x-1][y-1]-1]-1]-1]-1][y-1])%n
                            keepgoing=True
                        if M[x-1][a-1][b-1]!=Z and M[b-1][S[x-1][D[a-1][b-1]-1]-1][y-1]!=Z and M[S[x-1][D[a-1][b-1]-1]-1][a-1][S[b-1][D[S[x-1][D[a-1][b-1]-1]-1][y-1]-1]-1]!=Z and M[b-1][x-1][y-1]!=Z and M[x-1][a-1][S[b-1][D[x-1][y-1]-1]-1]==Z and M[S[b-1][D[x-1][y-1]-1]-1][S[x-1][D[a-1][S[b-1][D[x-1][y-1]-1]-1]-1]-1][y-1]!=Z:
                            M[x-1][a-1][S[b-1][D[x-1][y-1]-1]-1]= (M[b-1][S[x-1][D[a-1][b-1]-1]-1][y-1] +M[x-1][a-1][b-1] +M[S[x-1][D[a-1][b-1]-1]-1][a-1][S[b-1][D[S[x-1][D[a-1][b-1]-1]-1][y-1]-1]-1] -M[b-1][x-1][y-1]-M[S[b-1][D[x-1][y-1]-1]-1][S[x-1][D[a-1][S[b-1][D[x-1][y-1]-1]-1]-1]-1][y-1])%n
                            keepgoing=True
                        if M[x-1][a-1][b-1]!=Z and M[b-1][S[x-1][D[a-1][b-1]-1]-1][y-1]!=Z and M[S[x-1][D[a-1][b-1]-1]-1][a-1][S[b-1][D[S[x-1][D[a-1][b-1]-1]-1][y-1]-1]-1]!=Z and M[b-1][x-1][y-1]!=Z and M[x-1][a-1][S[b-1][D[x-1][y-1]-1]-1]!=Z and M[S[b-1][D[x-1][y-1]-1]-1][S[x-1][D[a-1][S[b-1][D[x-1][y-1]-1]-1]-1]-1][y-1]==Z:
                            M[S[b-1][D[x-1][y-1]-1]-1][S[x-1][D[a-1][S[b-1][D[x-1][y-1]-1]-1]-1]-1][y-1]= (M[b-1][S[x-1][D[a-1][b-1]-1]-1][y-1] +M[x-1][a-1][b-1] +M[S[x-1][D[a-1][b-1]-1]-1][a-1][S[b-1][D[S[x-1][D[a-1][b-1]-1]-1][y-1]-1]-1] -M[x-1][a-1][S[b-1][D[x-1][y-1]-1]-1]-M[b-1][x-1][y-1])%n
                            keepgoing=True
                        if M[x-1][a-1][b-1]!=Z and M[b-1][S[x-1][D[a-1][b-1]-1]-1][y-1]!=Z and M[S[x-1][D[a-1][b-1]-1]-1][a-1][S[b-1][D[S[x-1][D[a-1][b-1]-1]-1][y-1]-1]-1]!=Z and M[b-1][x-1][y-1]!=Z and M[x-1][a-1][S[b-1][D[x-1][y-1]-1]-1]!=Z and M[S[b-1][D[x-1][y-1]-1]-1][S[x-1][D[a-1][S[b-1][D[x-1][y-1]-1]-1]-1]-1][y-1]!=Z:
                            if (M[S[b-1][D[x-1][y-1]-1]-1][S[x-1][D[a-1][S[b-1][D[x-1][y-1]-1]-1]-1]-1][y-1])%n != (M[b-1][S[x-1][D[a-1][b-1]-1]-1][y-1] +M[x-1][a-1][b-1] +M[S[x-1][D[a-1][b-1]-1]-1][a-1][S[b-1][D[S[x-1][D[a-1][b-1]-1]-1][y-1]-1]-1] -M[x-1][a-1][S[b-1][D[x-1][y-1]-1]-1]-M[b-1][x-1][y-1])%n:
                                return False
    return tm2(M)

def bwdisplay(T):
    ### display Boltzmann weight as sum of char. functions
    out=[]
    for x in range(1,len(T[0])+1):
        for y in range(1,len(T[0])+1):
            for z in range(1,len(T[0])+1):
                if T[x-1][y-1][z-1]!=0:
                    out.append([T[x-1][y-1][z-1],x,y,z])
    return out


def bqsbwinv(g,X,bw,n):
    ### compute biquasile boltzmann weight invariant
    u=Symbol('u')
    PD=gauss2pd(g)
    H=pdcolors(X,PD)
    out=0
    for f in H:
        out=out+u**BW(PD,f,bw,n)
    return out


def BW(PD,f,bw,n):
    ### compute boltzmann weight of diagram
    out=0
    L,R=f[0],f[1]
    for x in PD:
        if x[0]==1:
            out=(out+bw[R[x[1]-1]-1][L[x[1]-1]-1][R[x[4]-1]-1])%n
        if x[0]==-1:
            out=(out-bw[L[x[4]-1]-1][L[x[3]-1]-1][R[x[4]-1]-1])%n
    return out

def phitester(bqs_matrix,n): 
### calculates which biquasiles are worth looking at
###
### args: 
### bqs = the entire matrix of biquasiles
### n = order of ring for bw 
	x = 0; 
	while x <= len(bqs_matrix):
		length = len(bwlist(bqs_matrix[x],n))
		if length !=1: 
			print("Biquasile index is " + str(x), "Length is " + str(length))
		else:
			print(str(x) + " biquasile indexes at 1")	
		x = x+1


def bwknotlist(bqs, bw, n):
	### calculates bqsbwinv for multiple knots at a time
	### 
	### args: 
	### bqs = an elt. of the breducelist
	### bw = an elt. of the bwlist
	### n = the prime for bw 
	L = []
	for x in gknotlist: 
		L.append([x,bqsbwinv(gknot[x],bqs,bw,n)])
	return L

def bwlinklist(bqs, bw, n):
	### calculates bqsbwinv for multiple knots at a time
	### 
	### args: 
	### bqs = an elt. of the breducelist
	### bw = an elt. of the bwlist
	### n = the prime for bw 
	L = []
	for x in glinklist: 
		L.append([x,bqsbwinv(glink[x],bqs,bw,n)])
	return L

def isEnhanceLink(bqs, bw, n): 
	### checks if enhnacement is trivial 
	### 
	### args: 
	### bqs = an elt. of the breducelist
	### bw = an elt. of the bwlist
	### n = the prime for bw 
	### returns:
	### link if enhanced by BW
	### nothing if enhnacement trivial 
	out = []
	for x in glinklist: 
		if type(bqsbwinv(glink[x],bqs,bw,n)) != Integer: 
			out.append("Biquasile is " + str(bqs) + " BW is " + str(bw) + " link is " + str(x))
	return out		

def isEnhanceKnot(bqs, bw, n): 
	### checks if enhancement is trivial
	### 
	### args: 
	### bqs = an elt. of the breducelist
	### bw = an elt. of the bwlist
	### n = the prime for bw 
	### returns:
	### knot if enhanced by BW
	### nothing if enhnacement trivial 
	out = []
	for x in gknotlist: 
		if type(bqsbwinv(gknot[x],bqs,bw,n)) != Integer: 
			out.append("Biquasile is " + str(bqs) + " BW is " + str(bw) + " knot is " + str(x))
	return out		

def bwknotrefrev(bqs, bw, n, g):
	### calculates bqsbwinv for the reverse and reflections of knots
	### 
	### args: 
	### bqs = an elt. of the breducelist
	### bw = an elt. of the bwlist
	### n = the prime for bw 
	### g = the gauss code for the knot
	L = []
	x = 0
	names = ["reg", "mirror", "reverse", "both"]
	while x < len(transformGC(g)): 
		L.append([names[x],bqsbwinv(transformGC(g)[x],bqs,bw,n)])
		x = x + 1 
	return L

def bwknotrefrevtester(bqs, bw, n, g):
	### calculates if any of the boltzmann weight differentiate between 
	### mirror images and reflections
	### 
	### returns:
	### 0 if there is no difference between all the knots 
	###
	### args: 
	### bqs = an elt. of the breducelist
	### bw = an elt. of the bwlist
	### n = the prime for bw 
	### g = the Gauss code for the knot 
	x = 0
	if not(bqsbwinv(transformGC(g)[0],bqs,bw,n) == bqsbwinv(transformGC(g)[1],bqs,bw,n) == bqsbwinv(transformGC(g)[2],bqs,bw,n) == bqsbwinv(transformGC(g)[3],bqs,bw,n)) :
		x = x+1
		print(g)
	return x

def bwknotloop(bqs, bw, n):
	### calculates if any of the boltzmann weight differentiate between 
	### mirror images and reflections
	### 
	### returns:
	### 0 if there is no difference between all the knots 
	###
	### args: 
	### bqs = an elt. of the breducelist
	### bw = an elt. of the bwlist
	### n = the prime for bw 
	z = 0
	for g in gknotlist: 
		z + bwknotrefrevtester(bqs, bw, n, g)
	return z

def bwlinkrefrev(bqs, bw, n, g):
	### calculates bqsbwinv for the reverse and reflections of knots
	### 
	### args: 
	### bqs = an elt. of the breducelist
	### bw = an elt. of the bwlist
	### n = the prime for bw 
	### g = the gauss code for the knot
	L = []
	x = 0
	names = ["reg", "mirror", "reverse", "both"]
	while x < len(transformGC2(g)): 
		L.append([names[x],bqsbwinv(transformGC2(g)[x],bqs,bw,n)])
		x = x + 1 
	return L

def bwknotlistrefrev(bqs, bw, n):
	### calculates bwknotrefrev for multiple knots at a time
	### 
	### args: 
	### bqs = an elt. of the breducelist
	### bw = an elt. of the bwlist
	### n = the prime for bw 
	L = []
	for x in gknotlist: 
		L.append([x,bwknotrefrev(bqs, bw, n, x)])
	return L

def bwlinklistrefrev(bqs, bw, n):
	### calculates bwknotrefrev for multiple knots at a time
	### 
	### args: 
	### bqs = an elt. of the breducelist
	### bw = an elt. of the bwlist
	### n = the prime for bw 
	L = []
	for x in glinklist: 
		L.append([x,bwlinkrefrev(bqs, bw, n, x)])
	return L




##################
# Alexander Stuff
##################

def zeromatrix(n):
    ### return n by n zero matrix### 
    out = []
    for i in range(0,n): 
        r = []
        for j in range(0,n): r.append(0)
        out.append(r)
    return out



def abqs(d,s,n,m):
    ### Alex Biquasile
    out=[zeromatrix(m),zeromatrix(m)]
    for x in range(1,m+1):
        for y in range(1,m+1):
            out[0][x-1][y-1]=(-d*s*n*n*x+n*y)%m
            out[1][x-1][y-1]=(d*x+s*y)%m
    for x in range(0,m):
        for y in range(0,m):
            for j in [0,1]:
                if out[j][x][y]==0:
                    out[j][x][y]=m
    return out

def abqslist(p,m):
	### Returns a list of all possible Alexander Biquasiles in Z_n
	### args:
	### p = p of Z/pZ
	### m = rank of matrix 

	out=[]
	d = 0
	while d < p: 
		s = 0
		while s < p: 
			n = 0 
			while n < p: 
				out.append(abqs(d,s,n,m))
				n = n+1
			s = s+1	
		d = d +1 
	return out	

def qhbw(X,w):
    ### homogeneous quadratic Boltzmann weight
    xx,xy,yy,xz,yz,zz=w[0],w[1],w[2],w[3],w[4],w[5]
    S,D,m=X[0],X[1],len(X[1])
    SRD,SLD,DRD,DLD=rdmatrix(S),ldmatrix(S),rdmatrix(D),ldmatrix(D)
    out=[]
    for j in range(1,len(S)+1):
        out.append(lm(zeromatrix(len(S))))
    for x in range(1,m+1):
        for y in range(1,m+1):
            for z in range(1,m+1):
                out[x-1][y-1][z-1]=(xx*x*x+xy*x*y+xz*x*z+yy*y*y+yz*y*z+zz*z*z)%m
    for x in range(1,len(S)+1):
        for a in range(1,len(S)+1):
            if out[x-1][a-1][ldiv(D,a,ldiv(S,x,x))-1]!=0: return False
            if out[x-1][rdiv(D,ldiv(S,x,x),a)-1][a-1]!=0: return False
    if not bwfill(X,m,out): return False
    return out


def qbw(X,w):
    ### quadratic Boltzmann weight
    xx,xy,yy,xz,yz,zz,A,B,C=w[0],w[1],w[2],w[3],w[4],w[5],w[6],w[7],w[8]
    S,D,m=X[0],X[1],len(X[1])
    SRD,SLD,DRD,DLD=rdmatrix(S),ldmatrix(S),rdmatrix(D),ldmatrix(D)
    out=[]
    for j in range(1,len(S)+1):
        out.append(lm(zeromatrix(len(S))))
    for x in range(1,m+1):
        for y in range(1,m+1):
            for z in range(1,m+1):
                out[x-1][y-1][z-1]=(xx*x*x+xy*x*y+xz*x*z+yy*y*y+yz*y*z+zz*z*z+A*x+B*y+C*z)%m
    for x in range(1,len(S)+1):
        for a in range(1,len(S)+1):
            if out[x-1][a-1][ldiv(D,a,ldiv(S,x,x))-1]!=0: return False
            if out[x-1][rdiv(D,ldiv(S,x,x),a)-1][a-1]!=0: return False
    if not bwfill(X,m,out): return False
    return out


def znm(n,m):
    ### generate (Zn)^m### 
    if m==1:
        out=[]
        for i in range(0,n):
            out.append([i])
        return out
    f = []
    for j in range(0,m):
        f.append(0)
    L = []
    while f:
       L.append(tuple(f))
       f = lexinc(f,n)
    return L
  

def lexinc(v,n):
    ### increment in dictionary order mod n### 
    p = len(v)-1   
    while p != -1 and v[p] == n-1:
        p = p-1
    if p == -1: 
        return False
    else:
        w = list(v)
        w[p] = w[p]+1 % n
        for j in range(p+1,len(w)):
            w[j] = 0
    return w



def qhbwlist(X):
    ### list quadratic Boltzmann weights
    out=[]
    for j in znm(len(X[0]),6):
        if qhbw(X,j):
            out.append(j)
    return out
    
def qbwlistN(X,N):
    ### list quadratic Boltzmann weights
    out=[]
    for j in znm(len(X[0]),9):
        if qbw(X,j):
            out.append(j)
            if len(out)>N: return out
    return out

def qbwlist(X):
    ### list quadratic Boltzmann weights
    out=[]
    for j in znm(len(X[0]),9):
        if qbw(X,j):
            out.append(j)
    return out

def lbwlist(X):
    ### list linear Boltzmann weights
    out=[]
    for j in znm(len(X[0]),3):
        if qbw(X,(0,0,0,0,0,0)+j):
            out.append(j)
    return out

def alexList(X, listfxn): 
	### returns all the non-trivial Alexander Biquasiles and their BW's
	### args: 
	### X = the list of matrices
	### lisfxn = the bw function applied to each matrix
	### returns:
	### out[0] = the list of biquasiles 
	### out[1] = the list of corresponding bw's
	j = 0
	bw =[]
	phi=[]
	out=[bw,phi]
	while j < len(X): 
		if len(listfxn(X[j])) != 1:
			bw.append(X[j])
			phi.append(listfxn(X[j]))

		j = j+1	

	return out


def bwify(X,phi):  		
	### Returns the bw values from the alpha/beta/gamma values
	### args:
	### X = biquasile matrix
	### phi = alpha, beta, gamma values for linear boltzmann
	return qbw(X, (0,0,0,0,0,0)+ phi)

def bwifyM(alex):
	### returns all the bw values from each Alexander biquasile 
	### args: 
	### alex = output of alexList 
	### returns:
	### out[0] = the list of biquasiles 
	### out[1] = the list of corresponding bw's

	bw = alex[0]
	phi = []
	out = [bw, phi]

	i = 0
	while i < len(alex[1]): 
		j = 0
		subcat = []
		while j < len(alex[1][i]): 
			subcat.append(bwify(alex[0][i], alex[1][i][j]))
			j = j+1
		phi.append(subcat)	
		i = i+1	
	return out 	

def alexEnhanceLink(alex,p,file):
	### Checks if the linear BW is an enhancement
	### args:
	### alex = output of bwifyM
	### p = p of Z/pZ
	### file = output file
	### writes on to txt file: 
	### biquasiles enhanced by BW
	### if none, writes failures
	i = 0
	out =[]
	while i < len(alex[0]):
		j = 0 
		while j < len(alex[1][i]):
			print("Hi! I'm still running! At: " + str(datetime.datetime.now()))
			x = isEnhanceLink(alex[0][i],alex[1][i][j],p)
			if (x): 
				out.append(x) 
				print("Success!" + str(x))
				file.write("Success!" + str(x))
				file.flush()
			else: 
				file.write("No luck for" + str(i) + " BW " + str(j) )
				file.flush()
			j = j +1
		i = i +1 
	if not out:
		print("For case " + str(p) + " all links fail")
		file.write("For case " + str(p) + " all links fail")	
		file.close()
	else: 
		print("For case " + str(p) + " we have" + str(out))
		file.write("For case " + str(p) + " we have" + str(out))
		file.close()
	print("link execution for " + str(p) + " finished") 	

    
def alexEnhanceKnot(alex,p,file):
	### Checks if the linear BW is an enhancement
	### args:
	### alex = output of bwifyM
	### p = p of Z/pZ
	### file = output file 
	### writes on to txt file: 
	### biquasiles enhanced by BW
	### if none, writes failures
	i = 0
	out =[]
	while i < len(alex[0]):
		j = 0 
		while j < len(alex[1][i]):
			print("Hi! I'm still running! At: " + str(datetime.datetime.now()))
			x = isEnhanceKnot(alex[0][i],alex[1][i][j],p)
			if (x): 
				out.append(x) 
				print("Success!" + str(x))
				file.write("Success!" + str(x))
				file.flush()
			else: 
				file.write("No luck for" + str(i) + " BW " + str(j) )
				file.flush()	
			j = j +1
		i = i +1 
	if not out:
		print("For case " + str(p) + " all knots fail")
		file.write ("For case " + str(p) + " all knots fail")	
		file.close()
	else: 
		print("For case " + str(p) + " we have" + str(out))
		file.write("For case " + str(p) + " we have" + str(out))
		file.close()
	print("knot execution for " + str(p) + " finished") 		


def easyKnotting(p, f):
	### Checks if linear BW enhancement over 
	### Alexander biquasile is trivial
	### args: 
	### p = characertistic of ring
	### f = file 
	alexEnhanceKnot(bwifyM(alexList(abqslist(p,p),lbwlist)),p,f)

def easyLinking(p, f):
	### Checks if linear BW enhancement over 
	### Alexander biquasile is trivial
	### args: 
	### p = characertistic of ring
	### f = file 
	alexEnhanceLink(bwifyM(alexList(abqslist(p,p),lbwlist)),p,f)
