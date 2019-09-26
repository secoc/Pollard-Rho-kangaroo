# python 2.7 !
import time
import random
import gmpy2
# windows https://code.google.com/archive/p/gmpy/downloads

modulo = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F
order  = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
Gx = 0x79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798
Gy = 0x483ada7726a3c4655da4fbfc0e1108a8fd17b448a68554199c47d08ffb10d4b8

class Point:
    def __init__(self, x=0, y=0):
        self.x = x
        self.y = y

PG = Point(Gx,Gy)
Z = Point(0,0) # zero-point, infinite in real x,y - plane

# return (g, x, y) a*x + b*y = gcd(x, y)
def egcd(a, b):
    if a == 0:
        return (b, 0, 1)
    else:
        g, x, y = egcd(b % a, a)
        return (g, y - (b // a) * x, x)

def rev(b, n = modulo):
    while b < 0:
        b += modulo
    g, x, _ = egcd(b, n)
    if g == 1:
        return x % n
        
def mul2(P, p = modulo):
    R = Point()
    c = 3*P.x*P.x*rev(2*P.y, p) % p
    R.x = (c*c-2*P.x) % p
    R.y = (c*(P.x - R.x)-P.y) % p
    return R

def add(P, Q, p = modulo):
    R = Point()
    dx = Q.x - P.x
    dy = Q.y - P.y    
    c = dy * gmpy2.invert(dx, p) % p     
    #c = dy * rev(dx, p) % p     
    R.x = (c*c - P.x - Q.x) % p
    R.y = (c*(P.x - R.x) - P.y) % p
    return R # 6 sub, 3 mul, 1 inv

def mulk(k, P = PG, p = modulo):
    if k == 0: return Z
    elif k == 1: return P
    elif (k % 2 == 0):
        return mulk(k/2, mul2(P, p), p)
    else:
        return add(P, mulk( (k-1)/2, mul2(P, p), p), p)

def X2Y(X, p = modulo):
    if p % 4 != 3:
        print 'prime must be 3 modulo 4'
        return 0
    X = (X**3+7)%p 
    pw = (p + 1) / 4 
    Y = 1
    for w in range(256):
        if (pw >> w) & 1 == 1:
            tmp = X
            for k in range(w):
                tmp = (tmp**2)%p
            Y *= tmp
            Y %= p
    return Y

def comparator():
    A, Ak, B, Bk = [], [], [], []
    with open('tame.txt') as f:
        for line in f:
            L = line.split()
            a = int(L[0],16)
            b = int(L[1])
            A.append(a)
            Ak.append(b)
    with open('wild.txt') as f:
        for line in f: 
            L = line.split()
            a = int(L[0],16)
            b = int(L[1])
            B.append(a)
            Bk.append(b)
    result = list(set(A) & set(B))
    if len(result) > 0: 
        sol_kt = A.index(result[0])
        sol_kw = B.index(result[0])
        print 'total time: %.2f sec' % (time.time()-starttime)
        d = Ak[sol_kt] - Bk[sol_kw]
        print 'SOLVED:', d
        file = open("results.txt",'a')
        file.write(('%d'%(Ak[sol_kt] - Bk[sol_kw])) + "\n")
        file.write("---------------\n")
        file.close()
        return True
    else:
        return False

def check(P, Pindex, DP_rarity, file2save):
    if P.x % (DP_rarity) == 0:
        file = open(file2save,'a')
        file.write(('%064x %d'%(P.x,Pindex)) + "\n")
        file.close()
        return comparator()
    else:
        return False
    
P = [PG]
for k in range(255): P.append(mul2(P[k]))    
print 'P-table prepared'    

def search():
    global solved
    DP_rarity = 1 << ((problem -  2*kangoo_power)/2 - 2)
    hop_modulo = ((problem-1) / 2) + kangoo_power 
    T, t, dt = [], [], []
    W, w, dw = [], [], []
    for k in range(Nt):
        t.append((3 << (problem - 2)) + random.randint(1, (1 << (problem - 1))))#-(1 << (problem - 2)) )
        T.append(mulk(t[k]))
        dt.append(0)
    for k in range(Nw):
        w.append(random.randint(1, (1 << (problem - 1))))
        W.append(add(W0,mulk(w[k])))
        dw.append(0)
    print 'tame and wild herds are prepared'
    oldtime = time.time()
    starttime = oldtime
    Hops, Hops_old = 0, 0
    t0 = time.time()
    oldtime = time.time()
    starttime = oldtime
    while (1):
        for k in range(Nt):
            Hops += 1
            pw = T[k].x % hop_modulo
            dt[k] = 1 << pw
            solved = check(T[k], t[k], DP_rarity, "tame.txt")
            if solved: break
            t[k] += dt[k]
            T[k] = add(P[pw], T[k])
        if solved: break            
        for k in range(Nw):
            Hops += 1
            pw = W[k].x % hop_modulo
            dw[k] = 1 << pw
            solved = check(W[k], w[k], DP_rarity, "wild.txt")
            if solved: break
            w[k] += dw[k]
            W[k] = add(P[pw], W[k])
        if solved: break
        t1 = time.time()
        if (t1-t0) > 5:
            print '%.3f h/s'%((Hops-Hops_old)/(t1-t0))
            t0 = t1
            Hops_old = Hops
    hops_list.append(Hops)        
    print 'Hops:', Hops        
    return 'sol. time: %.2f sec' % (time.time()-starttime)    

problems = [\
    ('029d8c5d35231d75eb87fd2c5f05f65281ed9573dc41853288c62ee94eb2590b7a',16),\
    ('036ea839d22847ee1dce3bfc5b11f6cf785b0682db58c35b63d1342eb221c3490c',24),\
    ('0209c58240e50e3ba3f833c82655e8725c037a2294e14cf5d73a5df8d56159de69',32),\
    ('03a2efa402fd5268400c77c20e574ba86409ededee7c4020e4b9f0edbee53de0d4',40),\
    ('025e466e97ed0e7910d3d90ceb0332df48ddf67d456b9e7303b50a3d89de357336',44),\
    ('026ecabd2d22fdb737be21975ce9a694e108eb94f3649c586cc7461c8abf5da71a',45),\
    ('03f46f41027bbf44fafd6b059091b900dad41e6845b2241dc3254c7cdd3c5a16c6',50),\
    ('0230210c23b1a047bc9bdbb13448e67deddc108946de6de639bcc75d47c0216b1b',65),\
    ('03bcf7ce887ffca5e62c9cabbdb7ffa71dc183c52c04ff4ee5ee82e0c55c39d77b',105)]

problem = 32
for elem in problems:
    s, n = elem
    if problem == n: break
kangoo_power = 3
Nt = Nw = 2**kangoo_power
X = int(s, 16)
Y = X2Y(X % (2**256))
if Y % 2 != (X >> 256) % 2: Y = modulo - Y
X = X % (2**256)
W0 = Point(X,Y)
starttime = oldtime = time.time()
search_range = 2**(problem-1)
Hops = 0
random.seed()

hops_list = []
N_tests = 3

for k in range(N_tests):
    solved = False
    open("tame.txt",'w').close()
    open("wild.txt",'w').close()
    search()
M = sum(hops_list)*1.0 / len(hops_list)
D = sum((xi - M) ** 2 for xi in hops_list)*1.0 / len(hops_list)
print M, '+/-',  (D / (len(hops_list)-1))**0.5
print 'Average time to solve: %.2f sec' % ((time.time()-starttime)/N_tests)
