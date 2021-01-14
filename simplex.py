import numpy as np
from scipy import linalg

A = np.array([[1,1,-1,0,0],[4,3,0,-1,0],[2,3,0,0,1]])
c = np.array([1,2,0,0,0])
b = np.array([35,120,150])

unsolved = True
i = 1
cols = [1,2,3]
B = A[:, cols]
cj = c[cols]
      
      
while unsolved:
    col_txt = [x+1 for x in cols]
    n_cols = [y for y in range(0,5) if y not in cols]
    n_cols_txt = [y+1 for y in n_cols]
    print('Starting iteration k = ', i)
    print('Basis B = [A{},A{},A{}]\n'.format(col_txt[0],col_txt[1],col_txt[2]), B)
    Binv = linalg.inv(B)
    print('B inverse \n',Binv)
    xb = Binv@b
    print('\nSolution {}'.format(i))
    print('xb = [x{}, x{}, x{}] = '.format(col_txt[0],col_txt[1],col_txt[2]), xb)
    print('xn = [x{}, x{}] = [0,0]'.format(n_cols_txt[0], n_cols_txt[1]))
    cb = np.array([c[z] for z in cols])
    cn = [c[t] for t in n_cols]
    print('cb = [c{}, c{}, c{}] = '.format(col_txt[0],col_txt[1],col_txt[2]), cb)
    print('cn = [c{}, c{}] = '.format(n_cols_txt[0], n_cols_txt[1]), cn)
    print('\nCompute reduced costs')
    cjs = []
    dn = []
    enter_basis = True
    for j in n_cols:
        cj = c[j] - cb@Binv@A[:,j]
        print('c{} = '.format(j+1), cj)
        if cj < 0:
            print('x{} is a candidate to enter the basis'.format(j+1))
            if enter_basis:
                x_enter = j
                dn.append(1)
                enter_basis = False #Use blands rule to assign lowest index of negative costs
            else:
                dn.append(0)
        else:
            dn.append(0)
        cjs.append(cj)
    if (cj >0 ).all():
        print('Solution found as all directions will increase cost')
        unsolved = False
        continue
    print('x{} will enter the basis'.format(x_enter+1))
    print('\nCompute feasible direction d')
    
    db = -Binv@A[:,x_enter]
    print('Direction db:', db)
    print('Non basic direction dn:', dn)
    if all(d > 0 for d in db):
        print('Solution found as there are no reduced costs')
        unsolved = False
        continue
            
    print('\nMinimum Ratio Test')
    thetas = []
    
    for xi in range(len(xb)):
        if db[xi]<0:
            theta = -xb[xi]/db[xi]
            print('theta{0} = x{0}/-db{0} = '.format(cols[xi]+1), theta)
            thetas.append(theta)
    theta_star = min(thetas)
    exit_basis = cols[thetas.index(theta_star)]
    print('\nx{} will exit the basis'.format(exit_basis+1))
    
    new_cols = [col for col in cols if col != exit_basis]
    new_cols.append(x_enter)
    new_n_cols = [col for col in n_cols if col != x_enter]
    new_n_cols.append(exit_basis)
    
    cols = new_cols
    cols.sort()
    n_cols = new_n_cols
    n_cols.sort()
    print('\nRedoing iteration using basis columns', [ct+1 for ct in cols])
    print('New non basic columns:', [nt+1 for nt in n_cols])
    B = A[:, cols]
    i+=1
    print('=======================================')
    if i == 10: #Catch errors for infinite loop
        unsolved = False
