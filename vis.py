import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as anim
import os

###============================================

def cdview(inputfile):
    ID    = []
    GROUP = []
    X     = []
    Y     = []
    Z     = []
    with open(inputfile) as f:
        Line = [s.strip() for s in f.readlines()]
        xlimn = float(Line[0][8:])
        ylimn = float(Line[1][8:])
        xlimp = float(Line[2][8:])
        ylimp = float(Line[3][8:])
        zlimn = float(Line[4][8:])
        zlimp = float(Line[5][8:])
        for l in range(6,len(Line)):
            id = ''
            g  = ''
            x  = ''
            y  = ''
            z  = ''
            index = 0
            for s in Line[l]:
                if s == ' ':
                    index += 1
                    continue
                if index == 0:
                    id += s
                elif index == 1:
                    g += s
                elif index == 2:
                    x += s
                elif index == 3:
                    y += s
                elif index == 4:
                    z += s
            ID.append(int(id))
            GROUP.append(int(g))
            X.append(float(x))
            Y.append(float(y))
            Z.append(float(z))

    colors = ["g", "r", "c", "m", "y", "k", "darksalmon", "limegreen", "b"]

    num_groups = max(GROUP) + 1
    X2 = [[] for _ in range(num_groups)]
    Y2 = [[] for _ in range(num_groups)]
    for i in range(len(GROUP)):
        X2[GROUP[i]].append(X[i])
        Y2[GROUP[i]].append(Y[i])

    ax.cla()
    ax.set_xlim(xlimn,xlimp)
    ax.set_ylim(ylimn,ylimp)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_facecolor("black")
    ax.set_aspect('equal')
    for i in range(num_groups):
        X = np.array(X2[i])
        Y = np.array(Y2[i])
        img = ax.scatter(X, Y, c=colors[i], marker="8", s=7.0, alpha=0.95)
###=============================================

fig = plt.figure(facecolor='black')
ax = fig.add_subplot(111)
FRAMES = []
for filename in sorted(os.listdir(".")):
    if ".cdv" in filename:
        FRAMES.append(filename)
ani = anim.FuncAnimation(fig, cdview, frames=FRAMES, interval=50)
gifname = "cdview.gif"
print("exporting",gifname,"...")
# ani.save(gifname, writer="imagemagick")
ani.save(gifname, writer="pillow")

###=============================================