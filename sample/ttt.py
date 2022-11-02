# LAMMPSの3次元シミュレーション結果から、2次元の配置を抜き出す

def three_to_two(path, output):
    X = []
    Y = []
    Z = []
    VX = []
    VY = []
    VZ = []
    step = 0
    bounds = []
    with open(path) as f:
        Line = [s.strip() for s in f.readlines()]
        is_step   = False
        is_bounds = False
        is_atoms  = False
        is_init = False
        for l in Line:
            if ("ITEM: TIMESTEP" in l):
                is_step = True
                is_bounds = False
                is_atoms = False
                is_init = False
                continue
            elif ("ITEM: BOX BOUNDS" in l):
                is_step = False
                is_bounds = True
                is_atoms = False
                continue
            elif ("ITEM: ATOMS" in l):
                is_step = False
                is_bounds = False
                is_atoms = True
                continue

            if (is_step and l=="0"):
                is_init = True
                is_step = False
                continue
            
            if is_init:
                continue
            
            if is_step:
                step = l
                is_step = False
            elif is_bounds:
                x_min = ''
                x_max = ''
                index = 0
                for s in l:
                    if s==' ':
                        index += 1
                    elif index==0:
                        x_min += s
                    elif index==1:
                        x_max += s
                bounds.append([x_min, x_max])
            elif is_atoms:
                x = ''
                y = ''
                z = ''
                vx = ''
                vy = ''
                vz = ''
                index = 0
                for s in l:
                    if s == ' ':
                        index += 1
                        continue
                    if index == 0:
                        x += s
                    elif index == 1:
                        y += s
                    elif index == 2:
                        z += s
                    elif index == 3:
                        vx += s
                    elif index == 4:
                        vy += s
                    elif index == 5:
                        vz += s
                if (not (-0.3<float(z)<0.3)):
                    continue
                X.append(x)
                Y.append(y)
                Z.append(z)
                VX.append(vx)
                VY.append(vy)
                VZ.append(vz)

    with open(output, "w") as f:
        f.write("ITEM: TIMESTEP\n")
        f.write(str(step)+"\n")
        f.write("ITEM: NUMBER OF ATOMS\n")
        f.write(str(len(X))+"\n")
        f.write("ITEM: BOX BOUNDS pp pp pp\n")
        for e in bounds:
            f.write("{} {}\n".format(e[0],e[1]))
        f.write("ITEM: ATOMS x y z vx vy vz\n")
        for i in range(len(X)):
            f.write(X[i]+" ")
            f.write(Y[i]+" ")
            f.write(Z[i]+" ")
            f.write(VX[i]+" ")
            f.write(VY[i]+" ")
            f.write(VZ[i])
            f.write("\n")



three_to_two("droplet.dump", "drp2d.dump")
