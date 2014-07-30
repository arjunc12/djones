import matplotlib
matplotlib.use('TkAgg')
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation

data = pd.read_table('rs_h', header=None, sep=r"\s*", engine='python', names=['x', 'y', 'z'])

frames = np.array_split(data, 9)

frames2 = np.array([])

for i in range(len(frames)):
    f = frames[i]
    grouped = f.groupby('x')
    unX = grouped['x'].unique()
    #print np.concatenate(tuple(unX))
    f2 = pd.DataFrame(columns = np.concatenate(tuple(unX)))
    j = 1
    for name, group in grouped:
        f2[j] = np.array(group['z'])
        j += 1
    frames[i] = f2
'''    
for frame in frames:
    plt.pcolor(frame)
    plt.show()
'''

def main():
    numframes = 9
    numpoints = 75
    
    x, y, c = np.random.random((3, numpoints))

    fig = plt.figure()
    heat = plt.pcolor(pd.DataFrame())

    ani = animation.FuncAnimation(fig, update_plot, frames=xrange(numframes), 
                                  interval = 50, blit=True)
    #ani.save("movie.mp4")
    plt.show()

def update_plot(i):
    frame = frames[i]
    heat = plt.pcolor(frame)
    return heat,

main()