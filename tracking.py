__author__ = 'Stella'

import scipy.io as sio

def convertMatFile(filename='Cell0000625_track.mat'):
	a=sio.matlab.loadmat(filename)

	# convert nan's to None


def simulatedAnneal():

from random import random

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
	E=0
	return E

def cost(state):
	E=0
	return E


convertMatFile()