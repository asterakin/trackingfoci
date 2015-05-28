__author__ = ['Stella', 'phil0']

from scipy import misc, io
import numpy as np
from random import random
import matplotlib.pyplot as plt


MIN_SCORE = 3
ALLOW_SPLITS = False
ALLOW_MERGES = False
MAX_TIME_WINDOW = 5



# returns false if all nan's
def no_nans(biglist):
    for x in biglist:
        for y in x:
            if not np.isnan(y):
                return False
    return True

def convertMatFile(filename):
    celldata = io.matlab.loadmat(filename)
    global lifetime
    global elements
    # TODO: does celldata have the key 'xx' ? I couldn't get this function to work with any of the .mat files I tried.
    # only use simpleTrack.mat - the other cell files are the old ones that don't work yet.

    lifetime = len(celldata['xx'][0])
    elements = len(celldata['xx'])
    # TODO: Please don't use the underscore as a variable name. It has a speacial python meaning - but also you just shouldn't use it. (fixed)
    # i thought they used _ in variable names in python for some reason..
    spots = [list([[np.nan] for y in range(lifetime)]) for i in range(elements)]
    for spot in range(elements):
        for time in range(lifetime):
            if celldata['sc'][spot][time] > MIN_SCORE:
                spots[spot][time] = [celldata['xx'][spot][time],celldata['yy'][spot][time],celldata['sc'][spot][time]]

    # spots has the format Track1: [[x1,y1,score],[x2,y2,score],[x3,y4,score],... ]
    # Track 2:  [[],[],[]]
    print(spots)
    plot(spots)
    return spots

def run (filename = 'simpleTrack.mat'):
    spots = convertMatFile(filename)
    # create random initial state
    # just take the spots the way they are in the table initially and act like they are a tracks
    # TODO: is tracks supposed to be global? Do you mean deepcopy (copy and deepcopy must be imported)? (If not, why is this here?)
    # i was just trying to copy the spots to a new list so that i have a new
    # copy of the tracks that i will be modifying.
    tracks = spots.copy()

		tracks = initial_state(tracks)
		[final_state,cost] = sim_anneal(tracks)

def initial_state (tracks):
    # do nearest neighbour from both sides to find some initial tracks
    #pick first spot
    # TODO: Why is this here? This is assigning the name newtrack to the range function (probably not whatever you're trying to do).
    # Not done yet.

    #newtrack = range
    #for t in range(1,lifetime):
		return tracks; # for now just keep the random initial state
    #pass

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
        T = T * alpha
    return state, cost

# evaluates cost of state
# TODO: Please don't push unless the code actually compiles (fixed).
def cost(state):
    # go through each track
    distance_metric = [None] * elements
    for track in range(len(state)):
        # TODO: don't leave your debugging code in here when you push (fixed).
        for time in range(1,lifetime):
            if state[track][time][0] is not None and state[track][time-1][0] is not None:
                distance_metric = (state[track][time][0]**2 - state[track][time-1][0]**2)**(0.5) 
                + (state[track][time][1]**2 - state[track][time-1][1]**2)**(0.5)

    cost = distance_metric

    return cost

def neighbor(state):
    # TODO: why are we setting the new_state to 0? Are states supposed to be ints? Set it to None if it's an object placeholder.

    # make random change in one spots

    new_state=None
    return new_state


def acceptance_probability(old_cost,new_cost,T)
    pass

def plot(tracks):
    # TODO: shouldn't 'Cell0000625.png' be an argument we pass in to the function? Or is this your debug code?
    # still trying to figure out how we plot and display in python. this is just a test.
    cellpicture = misc.imread('Cell0000625.png')
    plt.imshow(cellpicture)
    # plot two tracks
    for track in range(len(tracks)):
        newplot = []
        #if not empty(tracks[track]):
        for x in range(lifetime):
            newplot.append(tracks[track][x][0])
        plt.plot(range(0,lifetime),newplot)
        plt.scatter(range(0,lifetime),newplot)

