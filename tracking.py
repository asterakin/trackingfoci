__author__ = ['Stella', 'phil0']

from scipy import misc, io
import numpy as np
from random import randint
import matplotlib
import matplotlib.pyplot as plt
import math
from random import random
from operator import itemgetter

MIN_SCORE = 3
ALLOW_SPLITS = False
ALLOW_MERGES = False
MAX_TIME_WINDOW = 3
MAX_JUMP = 5


# returns false if all nan's
def no_nans(biglist):
    for x in biglist:
        for y in x:
            if not np.isnan(y):
                return False
    return True


# returns spots of format
# spots has the format Track1: [[x1,y1,score,hashcode],[x2,y2,score,hashcode],[x3,y4,score,hashcode],... ]
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
                spots[spot][time] = [celldata['xx'][spot][time], celldata['yy'][spot][time], celldata['sc'][spot][time],str(time)+'_'+str(spot)]

    return spots


def run(filename='simpleTrack.mat'):
    spots = convertMatFile(filename)
    tracks = spots.copy()
    tracks = initial_state_2(tracks)
    plt.ion()
    plot(tracks)
    [final_state, cost] = sim_anneal(tracks)
    plot(final_state)
    print("THE END")



def distance (pointA, pointB):
    print(pointA)
    print(pointB)
    return ((pointA[0] - pointB[0])**2 + (pointA[1] - pointB[1])**2)**0.5




def initial_state3 (state):
    # take random points, make random tracks towards both sides
    # then try to connect random segments with simulated annealing
    return state;



def initial_state_2(state):
    return state

# TODO: START BY FINDING THE NEAREST NEIGHBOR
def initial_state(state):
    #starts = np.empty([elements,4])
    #starts[:] = np.NAN

    # from, to, distance
    distance_map = []

    itracks_starts = [list([[np.nan] for y in range(lifetime)]) for i in range(elements)]
    itracks_ends = [list([[np.nan] for y in range(lifetime)]) for i in range(elements)]

    #ends = np.empty([elements,4])
    #ends[:] = np.NAN

    # take starting points
    for track in range(elements):
        if np.isfinite(state[track][0][0]):
            itracks_starts[track][0] = state[track][0]
            itracks_ends[track][0] = state[track][lifetime-1]



    for time in range(0, lifetime-1):
        for track1 in range(elements):
            pointA = itracks_starts[track1][time]
            for track2 in range(elements):
                pointB = state[track2][time+1]
                if np.isfinite(pointA[0]) and np.isfinite(pointB[0]):
                    d = distance(pointA,pointB)
                    distance_map.append([track1,track2,d])


        # find min - not a good check should prob use a different one.
        while distance_map !=[]:
            minDistance = float("inf");
            for x in range(len(distance_map)):
                if distance_map[x][2] < minDistance:
                    mintrack1 = distance_map[x][0]
                    mintrack2 = distance_map[x][1]
                    minDistance = distance_map[x][2]

            itracks_starts[mintrack1][time+1] = state[mintrack2][time+1]
            #print(itracks_starts)
        # remove elements for distance_map
            j=0
            while j < len(distance_map):
                if distance_map[j][0] == mintrack1 or distance_map[j][1] == mintrack2:
                    del distance_map[j]
                else:
                    j=j+1


    plot (itracks_starts)
    return state;  # for now just keep the random initial state


def sim_anneal(state):
    old_cost = cost(state)
    T = 1.0
    T_min = 0.00001
    alpha = 0.9
    while T > T_min:
        i = 1
        while i <= 10000:
            new_state = neighbor(state)
            new_cost = cost(state)
            ap = acceptance_probability(old_cost, new_cost, T)
            #print(str(old_cost) + ' vs ' + str(new_cost))
            if ap > random():
                print("Accepted")
                state = new_state
                old_cost = new_cost
                #plot(new_state)
            i += 1
            T = T * alpha

    return state, old_cost  # TODO: calculate a neighbor state

# TODO : Figure out the non-nan values and use only those
def neighbor(state):
    # make random change in random number of spots
    # swap random range

    how_many_spots = randint(0, 10)

    #for spots in range(how_many_spots):
    timepoint = randint(0, lifetime - 1)
    track1 = randint(0, 1)  # should be elements-1 but there are too many nan tracks
    track2 = randint(0, 1)
    temp = state[track1][timepoint:timepoint+how_many_spots]
    state[track1][timepoint:timepoint+how_many_spots] = state[track2][timepoint:timepoint+how_many_spots]
    state[track2][timepoint:timepoint+how_many_spots] = temp
    #plot(state)

    return state


def cost(state):

    distance_metric = [0 for i in range(elements)]
    for track in range(elements):
        # print("Track" + str(track))
        for time in range(0, lifetime):
            if time > 0 and np.isfinite(state[track][time][0]) and np.isfinite(state[track][time - 1][0]):
                distance_metric[track] = distance_metric[track] +(((state[track][time][0] - state[track][time - 1][0]) ** 2) + ((state[track][time][1] - state[track][time - 1][1]) ** 2)) ** 0.5

    icost = 0
    for i in distance_metric:
        icost = icost + i
    return icost


# TODO: figure out how to calculate acceptance probability
def acceptance_probability(old_cost, new_cost, T):
    ap = math.exp(old_cost - new_cost) / T
    return ap


def plot(tracks):
    # TODO: shouldn't 'Cell0000625.png' be an argument we pass in to the function? Or is this your debug code?

    #fig = plt.figure()
    plt.clf()
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
