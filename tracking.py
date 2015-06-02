__author__ = ['Stella', 'phil0']

from scipy import misc, io
import numpy as np
from random import randint
import matplotlib
import matplotlib.pyplot as plt
#from matplotlib.pyplot import plot, draw, show
import math
from random import random
import random
from operator import itemgetter
from copy import deepcopy

MIN_SCORE = 3
ALLOW_SPLITS = True
BIRTH_PENALTY = 20
ALLOW_MERGES = False
MAX_TIME_WINDOW = 5
MAX_JUMP = 5
MAX_JUMP_FUNC = lambda step: step * MAX_JUMP / 2.0


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


# spots has the format Track1: [[x1,y1,score,hashcode],[x2,y2,score,hashcode],[x3,y4,score,hashcode],... ]
#  Track 2:  [[],[],[]]
def convertMatFile(filename):
    celldata = io.matlab.loadmat(filename)
    global lifetime
    global elements
    global maxl

    splits = [np.nan for y in range(6)]
    merges = [np.nan for y in range(6)]

    lifetime = len(celldata['xx'][0])
    elements = len(celldata['xx'])

    maxl = celldata['lx'][0][lifetime-1]

    state = [list([[np.nan] for y in range(lifetime)]) for i in range(elements)]
    for spot in range(elements):
        for time in range(lifetime):
            if celldata['sc'][spot][time] > MIN_SCORE:
                state[spot][time] = [celldata['xx'][spot][time], celldata['yy'][spot][time], celldata['sc'][spot][time]]

    return state,splits,merges


def run(filename):
    [state,splits,merges] = convertMatFile(filename)
    [state,splits,merges] = deepcopy([state,splits,merges])
    plot(state)
    [final_state, c] = sim_anneal(state,splits,merges)
    plot(final_state)
    print('done')


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


def sim_anneal(state,splits,merges):
    old_cost = cost(state,splits,merges)
    T = 10.0
    T_min = 0.1
    alpha = 0.97
    old_cost_plot = []
    new_cost_plot = []
    while T > T_min:
        i = 1
        while i <= 100:
            new_state,new_splits,new_merges = neighbor_switch_jumps(state,splits,merges)
            new_cost = cost(new_state,new_splits,new_merges)
            ap = acceptance_probability(old_cost, new_cost, T)
            print(str(new_cost) +'vs'+ str(old_cost))
            if ap > random.random():
                print('accepted')
                state = new_state
                splits=new_splits
                merges=new_merges
                plot(new_state)
                old_cost = new_cost
            i += 1
            T = T * alpha
    return state,splits,merges,old_cost



def neighbor(state,splits,merges):
    # take one of this random options for operators
    operatorlist = []
    startflag=False
    endflag=False

    if count_big_jumps(state)>0:
         operatorlist.append(1)

    for start in find_starts_ends(state)[0]:
        if start > 0:
            startflag=True
            if ALLOW_SPLITS :
                operatorlist.append(2)

    for end in find_starts_ends(state)[1]:
        if end < lifetime:
            endflag=True
            if ALLOW_MERGES :
                operatorlist.append(3)

    if startflag and endflag:
        operatorlist.append(4)

    operator = random.choice(operatorlist)

    # exchange edges
    if operator ==1:
        return neighbor_switch_jumps(state,splits,merges)

    elif operator ==2:
    # find start connect it to a middle - split
        return neighbor_split(state,splits,merges)
        print('split')
    elif operator ==3:
    # find end connect it to a middle - merge
        return neighbor_merge(state,splits,merges)
        print('merge')
    elif operator ==4:
    # find start connect it to an end - closing gap
        return neighbor_closegap(state,splits,merges)
        print('close')
    else:
        return neighbor_onespot(state,splits,merges)


def neighbor_merge(state,splits,merges):
    possible_merge = []
    ends = find_starts_ends(state)[1]

    for itrack in range(elements):
        if ends[itrack] < lifetime:
            possible_merges.append(itrack)

    mergetrack = random.choice(possible_merges)
    ends_split = find_first(state[mergetrack])

    parent_tracks = []

    for itrack in range(elements):
        if starts[itrack] < start_split:
            parent_tracks.append(itrack)

    parent=random.choice(parent_tracks)

    splits[splittrack] = parent

    return state,splits,merges





def neighbor_split(state,splits,merges):
    possible_splits = []

    starts = find_starts_ends(state)[0]

    for itrack in range(elements):
        if starts[itrack] > 0:
            possible_splits.append(itrack)

    splittrack = random.choice(possible_splits)
    print(splittrack)
    start_split = find_first(state[splittrack])

    parent_tracks = []

    for itrack in range(elements):
        if starts[itrack] < start_split:
            parent_tracks.append(itrack)

    parent=random.choice(parent_tracks)

    splits[splittrack] = parent

    return state,splits,merges



def neighbor_close_gap(state,splits,merges):
    return state

def neighbor_switch_jumps(state,splits,merges):
    big_jumps = find_big_jumps(state)
    goodTrks = good_tracks(state)
    track1=random.choice(goodTrks)
    goodTrks.remove(track1)
    track2 = random.choice(goodTrks)

    num_of_jumps = len(big_jumps[track1])
    if num_of_jumps >= 1:
        which_jump = randint(0,num_of_jumps-1)
        time_jump1 = big_jumps[track1][0]
        print('switching ' + str(track1) + 'with' + str(track2) + 'at time' + str(time_jump1))
        temp = state[track1][time_jump1:lifetime]
        state[track1][time_jump1:lifetime] = state[track2][time_jump1:lifetime]
        state[track2][time_jump1:lifetime] = temp

    return state,splits,merges




def neighbor_onespot(state):
    # make random change for one random spots
    timepoint = randint(0, lifetime - 2)
    goodTrks = good_tracks(state)
    track1=random.choice(goodTrks)
    goodTrks.remove(track1)
    track2=random.choice(goodTrks)
    temp = state[track1][timepoint]
    state[track1][timepoint] = state[track2][timepoint]
    state[track2][timepoint] = temp

    temp2 = state[track1][timepoint+1]
    state[track1][timepoint+1] = state[track2][timepoint+1]
    state[track2][timepoint+1] = temp2

    return state


def find_starts_ends(state):
    starts = [np.nan for i in range(elements)]
    ends = [np.nan for i in range(elements)]

    for track in range(elements):
        starts[track]=find_first(state[track])
        ends[track]=find_last(state[track])

    return starts,ends



def cost(state,splits,merges):
    distance_metric = [0 for i in range(elements)]
    for track in range(elements):
        for time in range(1, lifetime):
            if np.isfinite(state[track][time][0]) and np.isfinite(state[track][time - 1][0]):
                distance_metric[track] += euclidean_distance(state[track][time], state[track][time - 1])

    icost = 0
    for i in distance_metric:
        icost = icost + i

    big_jump_count =count_big_jumps(state)

    #split cost
    splitcost=0
    for sTrack in range(elements):
        if find_first(state[sTrack]) != 0 and np.isnan(splits[sTrack]): # penanlty for starts in the middle of the timeline
            splitcost += BIRTH_PENALTY
        elif np.isfinite(splits[sTrack]):
            start = find_first(state[sTrack])
            parent = splits[sTrack]
            splitcost += euclidean_distance(state[sTrack][start], state[parent][time - 1])

    icost +=splitcost


    icost = icost + big_jump_count*2

    return icost


def find_big_jumps(state):
    total_result = [[] for i in range(elements)]
    for track in range(elements):
        result=[]
        for time in range(1, lifetime):
            if not np.isnan(state[track][time-1][0]) and not np.isnan(state[track][time][0]):
                distance = euclidean_distance(state[track][time-1], state[track][time])
                if distance > MAX_JUMP:
                    result.append (time)
        total_result[track]=result
    return total_result


def count_big_jumps(state):
    count=0
    for track in state:
        for i in range(len(track)):
            if i + 1 < len(track) and not np.isnan(track[i][0]) and not np.isnan(track[i+1][0]):
                distance = euclidean_distance(track[i], track[i + 1])
                if distance > MAX_JUMP:
                    count += 1
    return count



# TODO: figure out how to calculate acceptance probability
def acceptance_probability(old_cost, new_cost, T):
    ap = math.exp((old_cost - new_cost) / T)
    return ap



def initial_plot (state):
    for track in range(len(state)):
        newplot = []
        # if not empty(tracks[track]):
        for x in range(lifetime):
            newplot.append(state[track][x][0])
        ax.scatter(range(0, lifetime), newplot)




def plot(state):
    plt.clf()
    plt.axis([0, lifetime, -(maxl/2), maxl/2])
    for track in range(len(state)):
        newplot = []
        # if not empty(tracks[track]):
        for x in range(lifetime):
            newplot.append(state[track][x][0])
        plt.plot(range(0, lifetime), newplot, '.-')
        plt.draw()
    plt.show()




# points are of the form [x, y, intensity]
def euclidean_distance(point1, point2):
    return pow(pow(point1[0] - point2[0], 2) + pow(point1[1] - point2[1], 2), 0.5)



[state,splits,merges] = convertMatFile('simpleTrack_Cell0000741.mat')
c=cost(state,splits,merges)
#run('simpleTrack_Cell0000741.mat')
[st,s,m]=neighbor_split(state,splits,merges)
newc=cost(st,s,m)
print(str(c) +  ' vs '+ str(newc))



[st,s,m]=neighbor_split(state,splits,merges)
newc=cost(st,s,m)
print(str(c) +  ' vs '+ str(newc))


[st,s,m]=neighbor_split(state,splits,merges)
newc=cost(st,s,m)
print(str(c) +  ' vs '+ str(newc))



[st,s,m]=neighbor_split(state,splits,merges)
newc=cost(st,s,m)
print(str(c) +  ' vs '+ str(newc))



