import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation

data = pd.read_table('rs_h', header=None, sep=r"\s*", engine='python')

frames = np.array_split(data, 9)

def main():
    numframes = 9
    numpoints = 75
    
    x, y, c = np.random.random((3, numpoints))

    fig = plt.figure()
    scat = plt.scatter(x, y, c=c)#, s=100)

    ani = animation.FuncAnimation(fig, update_plot, frames=xrange(numframes), 
                                  interval = 10)
    #ani.save("movie.mp4")
    plt.show()

def update_plot(i):
    frame = frames[i]
    scat = plt.scatter(frame[0], frame[1], c=frame[2])
    return scat,

main()
    