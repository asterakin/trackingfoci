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
from future import *

MIN_SCORE = 3
ALLOW_SPLITS = True
BIRTH_PENALTY = 10
DEATH_PENALTY = 10

ALLOW_MERGES = True
MAX_TIME_WINDOW = 5
MAX_JUMP = 4
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



    lifetime = len(celldata['xx'][0])
    elements = len(celldata['xx'])

    maxl = celldata['lx'][0][lifetime-1]

    splits = [np.nan for y in range(elements)]
    merges = [np.nan for y in range(elements)]

    state = [list([[np.nan] for y in range(lifetime)]) for i in range(elements)]
    for spot in range(elements):
        for time in range(lifetime):
            if celldata['sc'][spot][time] > MIN_SCORE:
                state[spot][time] = [celldata['xx'][spot][time], celldata['yy'][spot][time], celldata['sc'][spot][time]]

    return state,splits,merges


def run(filename):
    [initial_state,splits,merges] = convertMatFile(filename)
    state = deepcopy(initial_state)
    plt.ion()
    plot(state,splits,merges)
    [final_state,splits,merges,c] = sim_anneal(initial_state,splits,merges)
    plot(final_state,splits,merges)
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
    T = 1.0
    T_min = 0.01
    alpha = 0.97
    iterations = 300
    old_cost_plot = []
    new_cost_plot = []
    while T > T_min:
        i = 1
        while i <= iterations:
            if i < iterations/2 :
                [new_state,new_splits,new_merges] = neighbor_switch_jumps(state,splits,merges)
            else:
               new_state,new_splits,new_merges = neighbor_merge_split(state,splits,merges)
            new_cost = cost(new_state,new_splits,new_merges)
            ap = acceptance_probability(old_cost, new_cost, T)
            print('new cost: ' +str(new_cost) +'vs old cost: '+ str(old_cost))
            if ap > random.random():
                print('accepted')
                state = deepcopy(new_state)
                splits=new_splits
                merges=new_merges
                plot(new_state,splits,merges)
                old_cost = new_cost
            i += 1
            T = T * alpha

    return state,splits,merges,old_cost



def neighbor_merge_split(state,splits,merges):
    # take one of this random options for operators
    operatorlist = []

    for start in find_starts_ends(state)[0]:
        if start > 0 and ALLOW_SPLITS :
                operatorlist.append(1)

    for end in find_starts_ends(state)[1]:
        if end < lifetime and ALLOW_MERGES :
                operatorlist.append(2)

    operator = random.choice(operatorlist)

    # exchange edges
    if operator ==1:
    # find start connect it to a middle - split
        print('Trying a split')
        return neighbor_split(state,splits,merges)

    elif operator ==2:
    # find end connect it to a middle - merge
        print('Trying a merge')
        return neighbor_merge(state,splits,merges)


def neighbor_merge(state,splits,merges):
    mergetrack=[]
    possible_merges = []
    ends = find_starts_ends(state)[1]
    #print('Each track ends : ' + str(ends))

    for itrack in range(elements):
        if ends[itrack] < lifetime-3:
            possible_merges.append(itrack)

    if possible_merges!=[]:
        mergetrack = random.choice(possible_merges)
        end_time = find_last(state[mergetrack])

        parent_tracks = []
        print(mergetrack)
        for itrack in range(elements):
            if ends[itrack] > end_time and not np.isnan(state[itrack][end_time+1][0]):
                parent_tracks.append(itrack)

        if parent_tracks!=[]:
            parent=random.choice(parent_tracks)
            merges[mergetrack] = parent

    return state,splits,merges



def neighbor_split(state,splits,merges):
    possible_splits = []

    starts = find_starts_ends(state)[0]

    for itrack in range(elements):
        if starts[itrack] > 0:
            possible_splits.append(itrack)

    splittrack = random.choice(possible_splits)
    start_split = find_first(state[splittrack])

    parent_tracks = []

    for itrack in range(elements):
        if starts[itrack] < start_split and np.isfinite(state[itrack][start_split][0]):
            parent_tracks.append(itrack)

    if parent_tracks != []:
        parent=random.choice(parent_tracks)
        splits[splittrack] = parent

    return state,splits,merges


# finds lonely spots and removes them from the track
def neighbor_remove_spots(state,splits,merges):




    return state

def neighbor_switch_jumps(state,splits,merges):
    big_jumps = find_big_jumps(state)
    goodTrks = good_tracks(state)
    track1=random.choice(goodTrks)
    goodTrks.remove(track1)


    num_of_jumps = len(big_jumps[track1])
    if num_of_jumps >= 1:
        which_jump = randint(0,num_of_jumps-1)
        time_jump1 = big_jumps[track1][which_jump]

        track2=random.choice(goodTrks) # if one is not found put a random one?
        for next_track in goodTrks:
            for k in big_jumps[next_track]:
                if k == time_jump1:
                    track2=next_track
                    break

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

    icost = icost/10

    big_jump_count =count_big_jumps(state)


    splitcost=0
    for sTrack in range(elements):
        if find_first(state[sTrack]) != 0 and np.isnan(splits[sTrack]): # penanlty for starts in the middle of the timeline
            splitcost += BIRTH_PENALTY
        elif np.isfinite(splits[sTrack]):
            start = find_first(state[sTrack])
            parent = splits[sTrack]
            if np.isfinite(state[sTrack][start][0]) and np.isfinite(state[parent][time - 1][0]):
                splitcost += euclidean_distance(state[sTrack][start], state[parent][time - 1])
            else:
                splitcost += 100

    mergecost=0
    for sTrack in range(elements):
        if find_last(state[sTrack]) < lifetime-1 and np.isnan(merges[sTrack]): # penanlty for starts in the middle of the timeline
            mergecost += DEATH_PENALTY
        elif np.isfinite(merges[sTrack]):
            end = find_last(state[sTrack])
            parent = merges[sTrack]
            if end < lifetime -2 and np.isfinite(state[sTrack][end][0]) and np.isfinite(state[parent][end + 1][0]):
                mergecost += euclidean_distance(state[sTrack][end], state[parent][end + 1])
            else:
                mergecost += 100

    # find nan time gaps
    nancount=0
    for track in range(elements):
        start = find_first(state[track])
        end = find_last(state[track])
        if not np.isnan(start):
            for time in range(start,end):
                if np.isnan(state[track][time][0]):
                    nancount +=1



    icost +=splitcost
    icost +=mergecost
    icost+= big_jump_count*5

    icost = icost/10 +nancount
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



def acceptance_probability(old_cost, new_cost, T):
    ap = math.exp((old_cost - new_cost) / T)
    return ap






def plot(state,splits,merges):
    plt.clf()
    plt.axis([0, lifetime, -(maxl/2), maxl/2])
    for track in range(elements):
        newplot = []
        gaps=[]
        gapst=[]
        start = find_first(state[track])
        end = find_last(state[track])

        if np.isfinite(splits[track]):
            parent = splits[track]
            plt.plot([start-1,start],[state[parent][start-1][0],state[track][start][0]],'--')

        if np.isfinite(start) and np.isfinite(end):
            for t in range(start,end+1):
                newplot.append(state[track][t][0])
                if np.isnan(state[track][t][0]):
                    if np.isfinite(state[track][t-1][0]):
                        gapst.append(t-1)
                        gaps.append(state[track][t-1][0])
                    if np.isfinite(state[track][t+1][0]):
                        gapst.append(t+1)
                        gaps.append(state[track][t+1][0])

            plt.plot(range(start, end+1), newplot, '.-')
            plt.plot(gapst, gaps, '--')
            plt.draw()
    plt.show()

# points are of the form [x, y, intensity]
def euclidean_distance(point1, point2):
    return pow(pow(point1[0] - point2[0], 2) + pow(point1[1] - point2[1], 2), 0.5)

run('simpleTrack.mat')
