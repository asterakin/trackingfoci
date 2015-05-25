__author__ = ['Stella', 'phil0']

import scipy.io as sio
from scipy import misc
import numpy as np
from random import random
import matplotlib.pyplot as plt

# should return false if all nan's
# not working great yet, returns false if there is any nan
def empty (biglist):
    for x in biglist:
	for y in x:
	    if np.isnan(y):
		return False
	    else:
		return True
        '''
        return all(np.isnan(x) for x in biglist)
        try:
            return all( empty(x) for x in biglist )
        except TypeError:
            return False'''

def convertMatFile(filename):
    celldata =sio.matlab.loadmat(filename)
    global lifetime
    global elements
    lifetime = len(celldata['xx'][0])
    elements = len(celldata['xx'])
    spots = [list([[numpy.nan] for y in range(lifetime)]) for _ in range(elements)]
    for spot in range(elements):
        for time in range(lifetime):
            if celldata['sc'][spot][time] > min_score:
                spots[spot][time] = [celldata['xx'][spot][time],celldata['yy'][spot][time],celldata['sc'][spot][time]]
    return spots

def run (filename = 'Cell0000625_track.mat'):
    global min_score
    global allow_splits
    global allow_merges
    global max_time_window
    # set parameters
    min_score= 3
    allow_splits= 0
    allow_merges= 0
    max_time_window = 5
    spots = convertMatFile(filename)
    # create random initial state
    # just take the spots the way they are in the table initially and act like they are a track
    tracks=spots.copy()

def initial_state (tracks):
    # do nearest neighbour from both sides to find some initial tracks
    #pick first spot
    newtrack = range
    for t in range(1,lifetime):
        pass

def sim_anneal(state):
    old_cost = cost(state)
    T = 1.0
    T_min = 0.00001
    alpha = 0.9
    while T > T_min:
        i = 1
        while i <= 100:
            new_state = neighbor(state)
            new_cost = cost(state)
            ap = acceptance_probability(old_cost, new_cost, T)
            if ap > random():
                state = new_state
                old_cost = new_cost
            i += 1
        T = T*alpha
    return state, cost

# evaluates cost of state
def cost(state):
    # go through each track
    distance_metric = [None for y in range(elements)]
    for track in range(len(state)):
        print('track:' +str(track))
        for time in range(1,lifetime):
            print(time)
            if state[track][time][0] is not None and state[track][time-1][0] is not None:
                    distance_metric = (state[track][time][0]**2 - state[track][time-1][0]**2)**(0.5) 
                    + (state[track][time][1]**2 - state[track][time-1][1]**2)**(0.5)
    return None

def neighbor(state):
    new_state=0
    return new_state

def plot(state):
    cellpicture = misc.imread('Cell0000625.png')
    plt.imshow(cellpicture)
    # plot two tracks
    for track in range(len(tracks)):
        newplot = []
        if not empty(tracks[track]):
            for x in range(lifetime):
                newplot.append(tracks[track][x][0])
            plt.plot(range(0,lifetime),newplot)
            plt.scatter(range(0,lifetime),newplot)

#convertMatFile()
