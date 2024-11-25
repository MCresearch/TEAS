import numpy as np

choice = int(input('choose the input file (0:SPIN1_CHG, 1:SPIN2_CHG):'))
file = 'SPIN1_CHG' if choice==0 else 'SPIN2_CHG'

with open(file, 'r') as f:
    temp = f.readlines()
    a = float(temp[1].split()[0]) / 0.529177249
    cell = np.array([[float(_) for _ in temp[2].split()],
                     [float(_) for _ in temp[3].split()],
                     [float(_) for _ in temp[4].split()],])
    element = temp[5].split()
    natom = np.array([int(_) for _ in temp[6].split()])

    start = 8
    position = np.zeros((np.sum(natom), 3))
    for i in range(np.sum(natom)):
        for j in range(3):
            position[i,j] = float(temp[start + i].split()[j])
        position[i] = position[i,0] * cell[0] + position[i,1] * cell[1] + position[i,2] * cell[2]
    position *= a
    position[:,2] += (192-135)/ 192 * a * 5
    
    start += np.sum(natom)
    nspin = int(temp[start + 1])
    fermi = float(temp[start + 2].split()[0])
    fftdim = np.array([int(_) for _ in temp[start + 3].split()])
    nx = fftdim[0]
    ny = fftdim[1]
    nz = fftdim[2]
    N = fftdim[0] * fftdim[1] * fftdim[2] 
    
    nrow = N // 8
    if N % 8 != 0:
        nrow += 1
    print(N)
    print(nrow)

    den = np.zeros((nx,ny,nz))
    count = 0
    for i in range(start + 5, start + 5 + nrow):
        data = [float(_) for _ in temp[i].split()]
        for j in range(len(data)):
            ix = int(count % nx)
            iy = int(((count - ix)/nx) % ny)
            iz = int((((count - ix)/nx) - iy) / ny)
            # print(ix, iy, iz)
            den[ix,iy,iz] = data[j]
            count += 1

    
    with open(file+'-trans.cube', 'w') as cube:
        cube.write('Cubefile created from ABACUS SCF calculation\n')
        cube.write('{0} (nspin) {1} (fermi energy, in Ry)\n'.format(nspin, fermi))
        cube.write('{0} 0.0 0.0 0.0\n'.format(np.sum(natom)))
        cell = cell * a
        cell[0] /= nx
        cell[1] /= ny
        cell[2] /= nz
        cube.write('{0} {1:.6f} {2:.6f} {3:.6f}\n'.format(nx, *cell[0]))
        cube.write('{0} {1:.6f} {2:.6f} {3:.6f}\n'.format(ny, *cell[1]))
        cube.write('{0} {1:.6f} {2:.6f} {3:.6f}\n'.format(nz, *cell[2]))
        for eachline in position:
            cube.write('13 13 {0:.6f} {1:.6f} {2:.6f}\n'.format(*eachline))
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    cube.write(' %.11e' % den[i,j,k-(192-135)])
                    if k % 6 == 5 and k != nz - 1:
                        cube.write('\n')
                cube.write('\n')