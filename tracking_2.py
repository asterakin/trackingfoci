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
MAX_TIME_WINDOW = 5
MAX_JUMP = 10 
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
    distance_map = [] # from, to, distance
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
                    d = euclidean_distance(pointA,pointB)
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
    return itracks_starts

def sim_anneal(state):
    old_cost = cost(state)
    T = 0.001
    T_min = 0.00001
    alpha = 0.9
    old_cost_plot = []
    new_cost_plot = []
    while T > T_min:
        i = 1
        while i <= 500:
            print(i)
            new_state = neighbor_onespot(state)
            new_cost = cost(state)
            ap = acceptance_probability(old_cost, new_cost, T)
            print(str(new_cost) +'vs'+ str(old_cost))
            if ap > random.random():
                print('accepted')
                state = new_state
                plot(new_state)
                old_cost = new_cost
            i += 1
            T = T * alpha
    return state, old_cost



def neighbor_switch_jumps(state):
    big_jumps = find_big_jumps(state)
    goodTrks = good_tracks(state)
    track1=random.choice(goodTrks)
    track2 = random.choice(goodTrks)


    num_of_jumps = len(big_jumps[track1])
    if num_of_jumps >= 1:
        which_jump = randint(0,num_of_jumps-1)
        time_jump1 = big_jumps[track1][0]
        temp = state[track1][time_jump1:lifetime]
        state[track1][time_jump1:lifetime] = state[track2][time_jump1:lifetime]
        state[track2][time_jump1:lifetime] = temp



    '''elif num_of_jumps>2:
        print('num of jumps' + str(num_of_jumps))
        print('list of big jumps' + str(big_jumps))
        which_jump = randint(0,num_of_jumps-2)
        print('which jump : ' + str(which_jump))
        time_jump1 = big_jumps[track1][which_jump]
        time_jump2= big_jumps[track1][which_jump+1]

        temp = state[track1][time_jump1:time_jump2]
        state[track1][time_jump1:time_jump2] = state[track2][time_jump1:time_jump2]
        state[track2][time_jump1:time_jump2] = temp
    '''

    return state

def neighbor(state):
    # make random change in random number of spots
    # swap random range
    timepoint = randint(0, lifetime - 1)
    goodTrks = good_tracks(state)
    track1=random.choice(goodTrks)
    goodTrks.remove(track1)
    track2=random.choice(goodTrks)
    temp = state[track1][timepoint:lifetime]
    state[track1][timepoint:lifetime] = state[track2][timepoint:lifetime]
    state[track2][timepoint:lifetime] = temp

    return state


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
        for time in range(1, lifetime):
            if np.isfinite(state[track][time][0]) and np.isfinite(state[track][time - 1][0]):
                distance_metric[track] += euclidean_distance(state[track][time], state[track][time - 1])

    icost = 0
    for i in distance_metric:
        icost = icost + i


    return icost

def cost2(state):
    distance_metric = [0 for i in range(elements)]
    msd = [0 for i in range(elements)]
    distance_from_average = [0 for i in range(elements)]
    meanx=[0 for i in range(elements)]
    meany=[0 for i in range(elements)]
    countspots = [0 for i in range(elements)]
    for track in range(elements):
        for time in range(0, lifetime):
            if np.isfinite(state[track][time][0]):
                meanx[track] += state[track][time][0]
                meany[track] += state[track][time][1]
                countspots[track] +=1
            if time>0 and np.isfinite(state[track][time][0]) and np.isfinite(state[track][time - 1][0]):
                distance_metric[track] += euclidean_distance(state[track][time], state[track][time - 1])
                msd[track] += euclidean_distance(state[track][time], state[track][0])

    #get difference from mean position of each track
    for track in range(elements):
        meanpoint = [meanx[track],meany[track]]
        for time in range(0, lifetime):
            if np.isfinite(state[track][time][0]):
                distance_from_average[track] += euclidean_distance(state[track][time],meanpoint)

    big_cont =count_big_continuities(state)
    icost = big_cont[0] + distance_from_average[0]+distance_metric[0]+msd[0]
    icost = icost/1000

    '''icost = 0
    for i in distance_metric:
        icost = icost + i

    for y in distance_from_average:
        icost = icost + y

    big_cont =count_big_continuities(state)
    for x in big_cont:
        icost +=x

    for i in msd:
        icost = icost + i

    icost = icost / 1000'''

    return icost

def count_big_jumps(state):
    result = []
    for track in state:
        count = 0
        for i in range(len(track)):
            if i + 1 < len(track) and not np.isnan(track[i][0]) and not np.isnan(track[i+1][0]):
                distance = euclidean_distance(track[i], track[i + 1])
                if distance > MAX_JUMP:
                    count += 1
        result.append(count)
    return result

def count_big_continuities(state):
    result = []
    for track in state:
        total_count = 0
        count = 0
        for i in range(len(track)):
            if i + 1 < len(track) and not np.isnan(track[i][0]) and not np.isnan(track[i+1][0]):
                distance = euclidean_distance(track[i], track[i + 1])
                if distance > MAX_JUMP:
                    total_count += pow(count,2)
                    count = 0 
                else:
                    count += 1
        total_count += pow(count,2)
        result.append(total_count)
    return result

# TODO: figure out how to calculate acceptance probability
def acceptance_probability(old_cost, new_cost, T):
    ap = math.exp(old_cost - new_cost) / T
    return ap

def plot(tracks):
    plt.clf()
    for track in range(len(tracks)):
        newplot = []
        # if not empty(tracks[track]):
        for x in range(lifetime):
            newplot.append(tracks[track][x][0])
        plt.plot(range(0, lifetime), newplot, '.-')

    plt.show()

    #plt.plot(range(0, lifetime), newplot, marker='o', markersize=5, linewidth=3)
        
    '''
	#cellpicture = misc.imread('Cell0000625.png')
    #plt.imshow(cellpicture)
	# plot two tracks
    #plt.figure(1)
    #plt.clf()
    #plt.plot(range(0,lifetime),newplot)
    plt.scatter(range(0,lifetime),newplot)
	'''


def search(state, data, current_time): 
    pass

def phil_nn(state):
    result = [[i[0]]  + (len(state[0]) - 1) * [[np.nan]] for i in state] 
    current_time = 0
    while (current_time + 1 < len(state[0])):
        already_examined_point_a = []
        for track in range(len(result)):
            point_a = result[track][current_time]
            if not np.isnan(point_a[0]) and point_a not in already_examined_point_a:
                already_examined_point_a.append(point_a)
                points_to_draw_line_to = find_nn(point_a, state, current_time, current_time + 1)
                result = attach_nn(track, current_time, current_time + 1, points_to_draw_line_to, result) 
                # couldn't find neighbor within max_jump, look more than one step ahead
                if not len(points_to_draw_line_to):
                    for step in range(2, MAX_TIME_WINDOW + 1):
                        if current_time + step < len(state):
                            points_to_draw_line_to = find_nn(point_a, state, current_time, current_time + step)
                            result = attach_nn(track, current_time, current_time + step, points_to_draw_line_to, result)
        current_time += 1
    return result

def find_nn(point, state, current_time, end_time):
    points_b = []
    potential_neighbors = [state[j][end_time] for j in range(len(state)) if not np.isnan(state[j][end_time][0])] 
    for neighbor_point in potential_neighbors:
        this_distance = euclidean_distance(point, neighbor_point) 
        if this_distance < MAX_JUMP:
            points_b.append(neighbor_point)
    return points_b

def attach_nn(track, current_time, end_time, points_b, result):
    for point_b in points_b:
        if not np.isnan(result[track][end_time][0]):
            for track2 in range(len(result)):
                if np.isnan(result[track2][current_time][0]):
                    result[track2][current_time] = result[track][current_time]
                    result[track2][end_time] = point_b
        else:
            result[track][end_time] = point_b
    return result

# points are of the form [x, y, intensity]
def euclidean_distance(point1, point2):
    return pow(pow(point1[0] - point2[0], 2) + pow(point1[1] - point2[1], 2), 0.5)




run()
state=convertMatFile('simpleTrack_Cell0000942.mat')
plot(state)
print('hi')



def sim_anneal(state,splits,merges):
    old_cost = cost(state,splits,merges)
    T = 10.0
    T_min = 0.01
    alpha = 0.97
    iterations = 500
    old_cost_plot = []
    new_cost_plot = []
    oplist = [1,2]

    while T > T_min:
        i = 1
        if i <= iterations:
            if i < 4*iterations/5.0 :
                operator = choice(oplist)
                if operator==1:
                    [new_state,new_splits,new_merges] = neighbor_switch_jumps(state,splits,merges)
                else:
                    [new_state,new_splits,new_merges] = look_ahead(state,splits,merges)
            else:
                new_state,new_splits,new_merges = neighbor_merge_split(state,splits,merges)
            new_cost = cost(new_state,new_splits,new_merges)
            ap = acceptance_probability(old_cost, new_cost, T)
            print('new cost: ' +str(new_cost) +' vs old cost: '+ str(old_cost))
            if ap > random.random() and old_cost != new_cost:
                print('accepted')
                state = deepcopy(new_state)
                splits=new_splits
                merges=new_merges
                plot(state,splits,merges)
                old_cost = new_cost
            i += 1
            T = T * alpha
    return state,splits,merges,old_cost