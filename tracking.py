__author__ = 'Stella'

import scipy.io as sio
from scipy import misc
from random import random


def is_empty (biglist):
	try:
		return all( is_empty(x) for x in biglist )
	except TypeError:
		return False

def convertMatFile(filename):
	celldata =sio.matlab.loadmat(filename)
	global lifetime
	global elements
	lifetime = len(celldata['xx'][0])
	elements = len(celldata['xx'])

	spots = [list([[None] for y in range(lifetime)]) for _ in range(elements)]

	for spot in range(elements):
		for time in range(lifetime):
			if celldata['sc'][spot][time] > min_score:
				spots[spot][time] = [celldata['xx'][spot][time],celldata['yy'][spot][time],celldata['sc'][spot][time]]


  return spots


def run (filename = 'Cell0000625_track.mat'):
		global min_score
		global allowSplits
		global allowMerges
		# set parameters

		min_score= 3
		allowSplits= 0
		allowMerges= 0

    spots = convertMatFile(filename)

    # create random initial state
		# just take the spots the way they are in the table initially and act like they are a track
    tracks=spots.copy()


def initial_state (tracks)
	# do nearest neighbour from both sides to find some initial tracks







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
					distance_metric = (state[track][time][0]**2 - state[track][time-1][0]**2)**(0.5) + (state[track][time][1]**2 - state[track][time-1][1]**2)**(0.5)



			E =

	return E

def neighbor(state):
	new_state=0
	return new_state

import matplotlib.pyplot as plt
from scipy import misc

def plot(state)
	cellpicture = misc.imread('Cell0000625.png')
	plt.imshow(cellpicture)


convertMatFile()