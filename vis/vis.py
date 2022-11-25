import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as anim
import os
import re

###============================================
# Muted by P.Tol
colors = ["#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#CC6677", "#882255", "#AA4499"]
# alternative
# colors = ["g", "r", "c", "m", "y", "k", "darksalmon", "limegreen", "b"]



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



def plot_energy(inputfile):
    STEP      = []
    KINETIC   = []
    POTENTIAL = []
    TOTAL     = []
    with open(inputfile) as f:
        Line = [s.strip() for s in f.readlines()]
        for l in range(0,len(Line)):
            step      = ''
            kinetic   = ''
            potential = ''
            total     = ''
            index = 0
            for s in Line[l]:
                if s == ' ':
                    index += 1
                    continue
                if index == 0:
                    step += s
                elif index == 1:
                    kinetic += s
                elif index == 2:
                    potential += s
                elif index == 3:
                    total += s
            STEP.append(int(step))
            KINETIC.append(float(kinetic))
            POTENTIAL.append(float(potential))
            TOTAL.append(float(total))
    
    print("exporting energy.png ...")
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.cla()
    ax.plot(STEP, KINETIC, label='kinetic', c=colors[4])
    ax.plot(STEP, POTENTIAL, label='potential', c=colors[2])
    ax.plot(STEP, TOTAL, label='total', c=colors[0])
    ax.set_xlabel('step')
    ax.set_ylabel('energy')
    ax.legend()
    plt.savefig('energy.png')

    

def plot_time(inputfile, outputfile):
    STEP = []
    T1   = []
    T2   = []
    T3   = []
    with open(inputfile) as f:
        Line = [s.strip() for s in f.readlines()]
        for l in range(0,len(Line)):
            step = ''
            t1   = ''
            t2   = ''
            t3   = ''
            index = 0
            for s in Line[l]:
                if s == ' ':
                    index += 1
                    continue
                if index == 0:
                    step += s
                elif index == 1:
                    t1 += s
                elif index == 2:
                    t2 += s
                elif index == 3:
                    t3 += s
            STEP.append(int(step))
            T1.append(float(t1))
            T2.append(float(t2))
            T3.append(float(t3))
    new_STEP = []
    T4 = []
    T5 = []
    T6 = []
    t4 = 0
    t5 = 0
    t6 = 0
    interval = 5
    for i in STEP:
        t4 += T1[i-1]
        t5 += T2[i-1]
        t6 += T3[i-1]
        if (i%interval==0):
            new_STEP.append(i)
            T4.append(t4/interval)
            T5.append(t5/interval)
            T6.append(t6/interval)
            t4 = 0
            t5 = 0
            t6 = 0

    
    print("exporting", outputfile, "...")
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.cla()
    ax.scatter(new_STEP, T4, label='min.', c=colors[4], marker='^', alpha=0.5, s=16.0)
    ax.scatter(new_STEP, T5, label='max.', c=colors[2], marker='o', alpha=0.8, s=16.0)
    ax.scatter(new_STEP, T6, label='avg.', c=colors[0], marker='x', alpha=0.7, s=16.0)
    ax.set_xlabel('step')
    ax.set_ylabel('time [ms]')
    ax.legend()
    plt.savefig(outputfile)

###=============================================
path = os.getcwd()
if re.search("myMD.vis", path):
    pass
else:
    os.chdir('vis')

cdv_path = ".."
for filename in os.listdir(".."):
    if filename == "cdv":
        cdv_path = "../cdv"
        break



# energy plot
plot_energy("../energy.dat")
plot_time("../time_whole.dat", "time_whole.png")
plot_time("../time_net.dat", "time_net.png")
plot_time("../time_gross.dat", "time_gross.png")
plot_time("../time_sdd.dat", "time_sdd.png")


# .cdv animation
plt.close()
fig = plt.figure(facecolor='black')
ax = fig.add_subplot(111)
FRAMES = []
for filename in sorted(os.listdir(cdv_path)):
    if ".cdv" in filename:
        FRAMES.append("{}/{}".format(cdv_path, filename))
ani = anim.FuncAnimation(fig, cdview, frames=FRAMES, interval=50)
gifname = "cdview.gif"
print("exporting",gifname,"...")
ani.save(gifname, writer="imagemagick")
# ani.save(gifname, writer="pillow")

###=============================================
