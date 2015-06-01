__author__ = ['Stella', 'phil0']

from scipy import misc, io
import numpy as np
from random import randint
import matplotlib
import matplotlib.pyplot as plt
import math
from random import random
import random
from operator import itemgetter
from copy import deepcopy

MIN_SCORE = 3
ALLOW_SPLITS = False
ALLOW_MERGES = False
MAX_TIME_WINDOW = 3
MAX_JUMP = 5


# returns true if all nan's
def has_nans(biglist):
    for x in biglist:
        for y in x:
            if np.isfinite(y):
                return False
    return True

# returns the indexes of tracks that are not all empty/full of nanas
def good_tracks(state):
    goodtracks = []
    for track in range(elements):
        if not has_nans(state[track]):
            goodtracks.append(track)
    return goodtracks



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
                spots[spot][time] = [celldata['xx'][spot][time], celldata['yy'][spot][time], celldata['sc'][spot][time]]
    return spots


def run(filename='simpleTrack.mat'):
    spots = convertMatFile(filename)
    tracks = deepcopy(spots)
    #tracks = initial_state(tracks)



    plt.ion()
    plot(tracks)

    print(find_starts_ends(tracks))
    [final_state, c] = sim_anneal(tracks)
    plot(final_state)
    print('done')



def distance (pointA, pointB):
    return ((pointA[0] - pointB[0])**2 + (pointA[1] - pointB[1])**2)**0.5

def find_first (track):
    for i in range(len(track)):
        if np.isfinite(track[i][0]):
            return (i)
    return np.nan

def find_last (track):
    for i in reversed(range(len(track))):
        if np.isfinite(track[i][0]):
            return (i)
    return(np.nan)

def initial_state2 (state):
    # take random points, make random tracks towards both sides
    # then try to connect random segments with simulated annealing
    return state


# TODO: START BY FINDING THE NEAREST NEIGHBOR
def initial_state(state):

    # from, to, distance
    distance_map = []

    itracks_starts = [list([[np.nan] for y in range(lifetime)]) for i in range(elements)]
    itracks_ends = [list([[np.nan] for y in range(lifetime)]) for i in range(elements)]

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

            j=0
            while j < len(distance_map):
                if distance_map[j][0] == mintrack1 or distance_map[j][1] == mintrack2:
                    del distance_map[j]
                else:
                    j=j+1

    plot (itracks_starts)
    return state


def sim_anneal(state):
    old_cost = cost(state)
    T = 1.0
    T_min = 0.00001
    alpha = 0.9
    old_cost_plot = []
    new_cost_plot = []
    while T > T_min:
        i = 1
        while i <= 1000:
            new_state = neighbor(state)
            new_cost = cost(state)
            ap = acceptance_probability(old_cost, new_cost, T)
            print(str(new_cost) +'vs'+ str(old_cost))
            if ap > random():
                print('accepted')
                state = new_state
                plot(new_state)
                old_cost = new_cost
            i += 1
            T = T * alpha

    return state, old_cost  # TODO: calculate a neighbor state

# TODO : Figure out the non-nan values and use only those
def neighbor(state):
    # make random change in random number of spots
    # swap random range

    #for spots in range(how_many_spots):
    timepoint = randint(0, lifetime - 1)



    goodTrks = good_tracks(tracks)
    track1=random.choice(goodTrks)
    goodTrks.remove(track1)
    track2=random.choice(goodTrks)


    temp = state[track1][timepoint:lifetime]
    state[track1][timepoint:lifetime] = state[track2][timepoint:lifetime]
    state[track2][timepoint:lifetime] = temp
    #plot(state)

    return state


def neighbor_onespot(state):
    # make random change in random number of spots
    # swap random range

    #for spots in range(how_many_spots):
    timepoint = randint(0, lifetime - 1)
    track1 = 0 # should be elements-1 but there are too many nan tracks
    track2 = 1
    temp = state[track1][timepoint]
    state[track1][timepoint] = state[track2][timepoint]
    state[track2][timepoint] = temp
    #plot(state)

    return state



def find_starts_ends(state):
    starts = [np.nan for i in range(elements)]
    ends = [np.nan for i in range(elements)]

    for track in range(elements):
        starts[track]=find_first(state[track])
        ends[track]=find_last(state[track])

    return starts,ends


def make_random_connections (state):
    # pick random track
    [starts,ends] = find_starts_ends(state)


    track1 = randint(0, elements)  # should be elements-1 but there are too many nan tracks
    track2 = randint(0, elements)
    #if track1 != track2:
     #   start = find_first(state[track1]
        #if start > 0:
            # connect to track2 - check score now.
         #   print('hi')
    pass






def cost(state):


    distance_metric = [0 for i in range(elements)]
    for track in range(elements):
        for time in range(0, lifetime):
            if time > 0 and np.isfinite(state[track][time][0]) and np.isfinite(state[track][time - 1][0]):
                distance_metric[track] += euclidean_distance(state[track][time], state[track][time - 1])
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
    plt.clf()
    for track in range(len(tracks)):
        newplot = []
        # if not empty(tracks[track]):
        for x in range(lifetime):
            newplot.append(tracks[track][x][0])
        plt.plot(range(0, lifetime), newplot, '.-')

    plt.show()


def phil_nn(state):
    result = [[i[0]] for i in state if not np.isnan(i[0][0])] 
    current_time = 0
    while (current_time + 1 < len(state[0])):
        potential_neighbors = [state[j][current_time + 1] for j in range(len(state)) if not np.isnan(state[j][current_time + 1][0])] 
        for track in range(len(result)):
            min_distance = float('inf')
            min_point = None
            for neighbor_point in potential_neighbors:
                this_distance = euclidean_distance(result[track][current_time], neighbor_point) 
                if this_distance < min_distance:
                    min_distance = this_distance
                    min_point = neighbor_point
            result[track].append(min_point)
        current_time += 1
    return result

# points are of the form [x, y, intensity]
def euclidean_distance(point1, point2):
    return pow(pow(point1[0] - point2[0], 2) + pow(point1[1] - point2[1], 2), 0.5)


run()


