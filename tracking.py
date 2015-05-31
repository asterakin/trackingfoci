__author__ = ['Stella', 'phil0']

from scipy import misc, io
import numpy as np
from random import randint
import matplotlib
import matplotlib.pyplot as plt
import math
from random import random

MIN_SCORE = 3
ALLOW_SPLITS = False
ALLOW_MERGES = False
MAX_TIME_WINDOW = 3



# returns false if all nan's
def no_nans(biglist):
    for x in biglist:
        for y in x:
            if not np.isnan(y):
                return False
    return True


# returns spots of format
# spots has the format Track1: [[x1,y1,score],[x2,y2,score],[x3,y4,score],... ]
#  Track 2:  [[],[],[]]
def convertMatFile(filename):
    celldata = io.matlab.loadmat(filename)
    global lifetime
    global elements

    lifetime = len(celldata['xx'][0])
    elements = len(celldata['xx'])

    spots = [list([[np.nan] for y in range(lifetime)]) for i in range(elements)]
    for spot in range(elements):
        for time in range(lifetime):
            if celldata['sc'][spot][time] > MIN_SCORE:
                spots[spot][time] = [celldata['xx'][spot][time], celldata['yy'][spot][time], celldata['sc'][spot][time]]

    return spots


def run(filename='simpleTrack.mat'):
    spots = convertMatFile(filename)
    tracks = spots.copy()
    tracks = initial_state(tracks)
    plt.ion()
    plot(tracks)
    [final_state, cost] = sim_anneal(tracks)
    plot(final_state)
    print("THE END")


# TODO: START BY FINDING THE NEAREST NEIGHBOR
def initial_state(tracks):
    return tracks;  # for now just keep the random initial state


def sim_anneal(state):
    old_cost = cost(state)
    T = 10.0
    T_min = 0.00001
    alpha = 0.9
    while T > T_min:
        i = 1
        while i <= 100:
            new_state = neighbor(state)
            new_cost = cost(state)
            ap = acceptance_probability(old_cost, new_cost, T)
            print(str(old_cost) + " vs " + str(new_cost))

            if ap > random():
                print("Accepted")
                state = new_state
                old_cost = new_cost
                plot(new_state)
            i += 1
            T = T * alpha

    return state, cost  # TODO: calculate a neighbor state

# TODO : Figure out the non-nan values and use only those
def neighbor(state):
    # make random change in random number of spots
    how_many_spots = randint(0, 10)

    for spots in range(how_many_spots):
        timepoint = randint(0, lifetime - 1)
        track1 = randint(0, 1)  # should be elements-1 but there are too many nan tracks
        track2 = randint(0, 1)
        # exchange track1 with track2
        temp = state[track1][timepoint]
        state[track1][timepoint] = state[track2][timepoint]
        state[track2][timepoint] = temp

    return state


def cost(state):
    distance_metric = [0 for i in range(elements)]
    for track in range(elements):
        # print("Track" + str(track))
        for time in range(0, lifetime):
            if time > 0 and np.isfinite(state[track][time][0]) and np.isfinite(state[track][time - 1][0]):
                distance_metric[track] = distance_metric[track] +(((state[track][time][0] - state[track][time - 1][0]) ** 2) + ((state[track][time][1] - state[track][time - 1][1]) ** 2)) ** 0.5
            if time +1 < lifetime and np.isfinite(state[track][time][0]) and np.isfinite(state[track][time + 1][0]):
                distance_metric[track] = distance_metric[track] +(((state[track][time][0] - state[track][time + 1][0]) ** 2) + ((state[track][time][1] - state[track][time + 1][1]) ** 2)) ** 0.5

    cost = 0
    for i in distance_metric:
        cost = cost + i
    return cost


# TODO: figure out how to calculate acceptance probability
def acceptance_probability(old_cost, new_cost, T):
    ap = math.exp(old_cost - new_cost) / T
    return ap


def plot(tracks):
    # TODO: shouldn't 'Cell0000625.png' be an argument we pass in to the function? Or is this your debug code?

    fig = plt.figure()
    plt.close(fig)
    print("I am plotting")
    for track in range(len(tracks)):
        newplot = []
        # if not empty(tracks[track]):
        for x in range(lifetime):
            newplot.append(tracks[track][x][0])
        plt.plot(range(0, lifetime), newplot, '.-')

    plt.show()

    '''
		#cellpicture = misc.imread('Cell0000625.png')
    #plt.imshow(cellpicture)
		# plot two tracks
    #plt.figure(1)
    #plt.clf()
    #plt.plot(range(0,lifetime),newplot)
    plt.scatter(range(0,lifetime),newplot)
		'''


run()
