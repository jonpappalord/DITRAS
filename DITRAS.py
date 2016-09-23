#!/usr/bin/env python
# coding: utf-8

from random import sample, uniform, random
import cPickle as pickle
from math import sqrt, sin, cos, pi, asin, pow, ceil
from scipy.stats import expon
from collections import defaultdict
import numpy as np
import sys
import os
from os import path
import time
import csv

__author__ = "Luca Pappalardo and Filippo Simini"
__contact__ = "luca.pappalardo@isti.cnr.it"
__license__ = "BSD"
__credits__ = "Luca Pappalardo and Filippo Simini, " \
              "Modelling individual routines and spatio-temporal trajectories in human mobility, " \
              "http://arxiv.org/abs/1607.05952, " \
              "2016"

def timeit(method):
    """
    Compute the execution time of a function or a method.
    """
    def timed(*argst, **kwt):
        ts = time.time()
        result = method(*argst, **kwt)
        te = time.time()

        print '\n\t\ttime: %2.2f min' % ((te-ts)/60.0)
        return result

    return timed

@timeit
def compute_od_matrix(spatial_tessellation, filename='od_matrix.pkl'):
    """
    Compute a weighted origin destination matrix where an element A_{ij} is the
    probability p_{ij} of moving between two locations in the spatial tessellation
    given as input.

    Parameters
    ----------
    spatial_tessellation: dict
        a dictionary of location identifiers to a dictionary containing latitude,
        longitude and density of the location

    filename: string
        the name of the file where to store the computed origin destination matrix
        default: 'od_matrix.pkl'

    Returns
    -------
    od_matrix: numpy array
        a bidimensional numpy array describing the weighted origin destination matrix
    """

    sys.stdout.write('[Computing origin-destination matrix]\n')
    sys.stdout.flush()

    n = len(spatial_tessellation)
    old_p = 0
    bar_length = 20
    od_matrix = np.zeros( (n, n) )
    count = 1
    for id_i in spatial_tessellation:
        lat_i, lon_i, d_i = spatial_tessellation[id_i]['lat'], spatial_tessellation[id_i]['lon'], \
                            spatial_tessellation[id_i]['relevance']
        for id_j in spatial_tessellation:
            if id_j != id_i:
                lat_j, lon_j, d_j = spatial_tessellation[id_j]['lat'], spatial_tessellation[id_j]['lon'], \
                                    spatial_tessellation[id_j]['relevance']
                p_ij = (d_i * d_j) / (earth_distance((lat_i, lon_i), (lat_j, lon_j)) ** 2)
                od_matrix[id_i, id_j] = p_ij

        # normalization by row
        od_matrix[id_i] /= np.sum(od_matrix[id_i])

        # progress bar update
        percentage = int(float(count * 100) / n)
        if percentage > old_p:
            hashes = '-' * int(round(percentage/5))
            spaces = ' ' * (bar_length - len(hashes))
            sys.stdout.write("\rExec: [{0}] {1}%".format(hashes + spaces, int(round(percentage))))
            sys.stdout.flush()
            old_p = percentage
        count += 1

    pickle.dump(od_matrix, open(filename, 'wb'))
    return od_matrix

def earth_distance(lat_lng1, lat_lng2):
    """
    Compute the distance (in km) along earth between two latitude and longitude pairs

    Parameters
    ----------
    lat_lng1: tuple
        the first latitude and longitude pair
    lat_lng2: tuple
        the second latitude and longitude pair

    Returns
    -------
    float
        the distance along earth in km
    """
    lat1, lng1 = [l*pi/180 for l in lat_lng1]
    lat2, lng2 = [l*pi/180 for l in lat_lng2]
    dlat, dlng = lat1-lat2, lng1-lng2
    ds = 2 * asin(sqrt(sin(dlat/2.0) ** 2 + cos(lat1) * cos(lat2) * sin(dlng/2.0) ** 2))
    return 6371.01 * ds  # spherical earth...

def weighted_random_selection(weights):
    """
    Choose an index from the list of weights according to the numbers in the list
    
    Parameters
    ----------
    weights: list
        a list of weights (e.g., probabilities)
        
    Returns
    -------
    index: int
        the index of element chosen from the list
    """
    totals = []
    running_total = 0

    for w in weights:
        running_total += w
        totals.append(running_total)

    rnd = random() * running_total
    for index, total in enumerate(totals):
        if rnd < total:
            return index

def load_spatial_tessellation(filename='location2info', delimiter=','):
    """

    FARE CON PANDAS

    Load into a dictionary the locations and corresponding information (latitude, longitude, relevance)
    
    Parameters
    ----------
    filename: str
        the filename where the location info is stored
        
    Returns
    -------
    dict
        the dictionary of locations
    """
    spatial_tessellation = {}
    f = csv.reader(open(filename), delimiter=delimiter)
    f.next()  # delete header
    i = 0
    for line in f:
        spatial_tessellation[i] = {'lat': float(line[0]),
                                   'lon': float(line[1]),
                                   'relevance': int(line[2])}
        i += 1
    return spatial_tessellation

class TrajectoryGenerator(object):
    """
    Superclass describing a trajectory generator. A trajectory generator must implement two methods:
        - start_simulation, which specifies how to assign physical locations to a mobility diary
        - choose_location, which specifies how to choose the next location during the simulation
    """

    def __init__(self, name):
        """
        Constructor
        
        Parameters
        -----------
        name: string
            the name of the trajectory generator
        
        """
        self.name = name
        
    def choose_location(self):
        """
        It specifies how to choose the next location on the given spatial tessellation
        
        """
        return ""
    
    def start_simulation(self, spatial_tessellation, mobility_diary, od_matrix):
        """
        Specifies how the model works
        
        Parameters
        ----------
        spatial_tessellation: dict
            a dictionary of locations identifiers to location infos

        mobility_diary: list
            the mobility diary on which to base the generation of trajectory
        """
        return []
        
class dEPR(TrajectoryGenerator):
    """
    The d-EPR mobility model as described in:

        - Pappalardo et al.,
          Returners and Explorers dichotomy in Human Mobility,
          Nature Communications 6:8166 doi: 10.1038/ncomms9166, 2015

        - Pappalardo et al.,
          Human Mobility Modelling: exploration and preferential return meet the gravity model,
          http://dx.doi.org/10.1016/j.procs.2016.04.188, 2016.

    """
    def __init__(self, rho, gamma):
        self.name = 'd-EPR'
        self.rho = rho
        self.gamma = gamma
    
    def __preferential_return(self):
        """
        Choose the location the individual returns to, according to the visitation frequency to the already visited locations

        Returns
        -------
        next_location: int
            the identifier of the next location
        """
        index = weighted_random_selection(self.location2visits.values())
        next_location = self.location2visits.keys()[index]
        return next_location

    def __preferential_exploration(self, current_location):
        """
        Choose the new location the individual explores, according to the probabilities in matrix M.

        Parameters
        ----------
        current_location: int
            the identifier of the current location of the individual

        Returns
        -------
        next_location: int
            the identifer of the new location to explore
        """
        return weighted_random_selection(self.od_matrix[current_location])

    def choose_location(self):

        # initialize variables
        S = len(self.location2visits.keys())  # number of already visited locations
        
        if S == 0:
            self.home = self.__preferential_exploration(self.home)
            return self.home

        ## choose a probability to return o explore
        p_new = uniform(0, 1)

        if p_new <= self.rho * pow(S, -self.gamma):  # choose to return or explore
            # PREFERENTIAL EXPLORATION
            current_location = self.trajectory[-1]  # the last visited location
            return self.__preferential_exploration(current_location)

        else:
            # PREFERENTIAL RETURN
            return self.__preferential_return()


    def start_simulation(self, spatial_tessellation, mobility_diary, od_matrix):
        """
        Start the simulation of human mobility based on the mobility diary given in input

        Parameters
        ----------
        mobility_diary: list
            the mobility diary on which to base the simulation

        Returns
        -------


        """
        ## initialization of parameters
        self.trajectory = []
        self.location2visits = defaultdict(int)
        self.od_matrix = od_matrix
        self.spatial_tessellation = spatial_tessellation
        self.home = sample(self.spatial_tessellation.keys(), 1)[0]
        
        i = 0
        while i < len(mobility_diary):

            if mobility_diary[i] == 0:  # the agent is at home
                self.trajectory.append(self.home)
                self.location2visits[self.home] += 1
                i += 1

            else:  # the agent is not at home
                next_location = self.choose_location()
                self.trajectory.append(next_location)
                self.location2visits[next_location] += 1
                j = i + 1

                while j < len(mobility_diary) and mobility_diary[i] == mobility_diary[j]:
                    self.trajectory.append(next_location)
                    self.location2visits[next_location] += 1
                    j += 1
                i = j 

        return self.trajectory


class DiaryGenerator(object):
    """
    Superclass describing a temporal model. A temporal model must implement a start_simulation method.
    """
    def __init__(self):
        self.name = ""
    
    def start_simulation(self, diary_length):
        return []
        
class RandomDiary(DiaryGenerator):
    
    def __init__(self):
        self.name = "random_diary"
        
    def start_simulation(self, diary_length):
        return [i for i in range(1, diary_length + 1)]

    
class WaitingTimeDiary(DiaryGenerator):
    def __init__(self, beta=0.8, tau=17):
        self.name = "waiting_time"
        self.beta = beta
        self.tau = tau # in hours
        
    def __get_waiting_time(self):
        """
        Extract a waiting time from a power law with exponential cut-off distribution.
        The parameters of the distribution are taken from the paper:
        C. Song et al., Modelling the scaling properties of human mobility, Nature Physics 6, 818-823 (2010).

        ---
        To simulate a power law with exponential cut-off x^(-alpha) * exp(-lambda * x), we can generate an exponentially
        distributed random number U and then accept or reject it with probability p or 1-p respectively (i.e. accept if U < p
        or reject if U > p, where U is a uniform [0, 1] random variable), where p = (x/x_min)^(-alpha) and x_min=1.

        http://www.santafe.edu/aaronc/powerlaws/
        ---

        :return: float
            a waiting time chosen from the waiting time distribution
        """
        
        x = expon.rvs(1.0/self.tau)
        while pow(x, -(1 + self.beta)) < uniform(0.0, 1.0):
            x = expon.rvs(1.0/self.tau)

        return x
        
    def start_simulation(self, diary_length):
        i, total_count = 0, 1
        D = []
        while i < diary_length:
            waiting_time = int(ceil(self.__get_waiting_time()))
            for j in range(0, waiting_time):
                D.append(total_count)
            total_count += 1
            i += waiting_time
            
        return D[:diary_length]
            
    
class MD(DiaryGenerator):
    """
    The MD diary generator.
    Reference:
        - Pappalardo et al., Modelling individual routines and spatio-temporal trajectories in human mobility,
        http://arxiv.org/abs/1607.05952, 2016
    """
    def __init__(self, filename='diary_generator_1hour.pkl'):
        self.name = "MD"
        self.markov_model = pickle.load(open(filename, 'rb'))
        
        
    def start_simulation(self, diary_length):
        """
        Start the simulation of the mobility diary given the diary_length

        Parameters
        ----------
        diary_length: int
            the length of the diary in the time unit of the MD Markov chain

        Returns
        -------
        D: list
            the mobility diary (a list of 0 and 1)
        """
        V, i = [], 0
    
        prev_state = (i, 1) ## it starts from the typical location at midnight
        ### WARNING: this has to be changed to a more proper choice
        
        V.append(prev_state)

        while i < diary_length:
            
            h = i % 24 ## the hour of the day

            ## select the next state in the Markov chain
            index = weighted_random_selection(self.markov_model[prev_state].values()) 
            next_state = self.markov_model[prev_state].keys()[index]
            V.append(next_state)

            j = next_state[0]
            if j > h:  # we are in the same day
                i += j - h
            else:      # we are in the next day
                i += 24 - h + j

            prev_state = next_state

        ### now we translate the temporal diary into the the mobility diary
        prev, D, other_count = V[0], [], 1
        D.append(0) #### WARNING: this has to be changed to a more proper choice
        
        for v in V[1:]: ## scan all the states obtained and create the synthetic time series
            h, s = v
            h_prev, s_prev = prev

            if s == 1:        ## if in that hour she visits home
                D.append(0)
                other_count = 1
            else:             ## if in that hour she does NOT visit home
                
                if h > h_prev: ### we are in the same day
                    j = h - h_prev
                else:          ### we are in the next day
                    j = 24 - h_prev + h

                for i in range(0, j):
                    D.append(other_count)
                other_count += 1

            prev = v

        return D[0: diary_length]

@timeit
def load_od_matrix(od_matrix):
    return pickle.load(open(od_matrix, 'rb'))

class DITRAS(object):
    """
    The DITRAS (DIary-based TRAjectory Simulator) model as described in:

        - L. Pappalardo and F. Simini, Modelling individual routines and spatio-temporal trajectories in human mobility,
          http://arxiv.org/abs/1607.05952
    """
    def __init__(self, n_agents=10000, length=168, diary_generator=MD(),
                 trajectory_generator=dEPR(rho=0.6, gamma=0.21), filename='ditras_trajs.csv'):
        """
        Constructor

        Parameters
        ----------
        n_agents: int
            the number of agents to simulate

        length: int
            the length of the diary (in time units depending on the diary generator)

        diary_generator: DiaryGenerator object
            the diary generator

        trajectory_generator: TrajectoryGenerator object
            the trajectory generator

        filename: string
            the name of the file where to write the produced synthetic trajectories

        """
        self.n_agents = n_agents
        self.length = length
        self.diary_generator = diary_generator
        self.trajectory_generator = trajectory_generator
        self.filename = filename

    @timeit
    def start_simulation(self, spatial_tessellation, od_matrix=None):
        """
        Start the simulation of DITRAS on the spatial tessellation given as input

        Parameters
        ----------
        spatial_tessellation: dict
            a dictionary of location identifiers to a dictionary of location infos

        od_matrix: None or string
            if None, the weighted origin destination matrix (used by the trajectory generator d-EPR)
             will be computed based on the spatial tessellation given as input.
            if it is a string, it indicates the filename of a pickle object where a precomputed
             origin destination matrix is stored (e.g., od_matrix.pkl).
        """
        sys.stdout.write('\n[Loading OD matrix]\n')
        self.od_matrix = load_od_matrix(od_matrix)

        bar_length = 20
        old_p = 0

        sys.stdout.write('\n[DITRAS simulation]\n')
        sys.stdout.write('\t%s agents\n\t%s time slots\n\tdiary gen: %s\n\ttraj gen: %s\n\n'
                         % (self.n_agents, self.length,
                           self.diary_generator.__class__.__name__, self.trajectory_generator.__class__.__name__))
        sys.stdout.flush()

        f = open(self.filename, 'w')
        f.write('user,location,time_slot\n')

        n_rows = 0
        count = 1
        for agent_id in range(0, self.n_agents):
            mobility_diary = self.diary_generator.start_simulation(self.length)

            synthetic_trajectory = \
                self.trajectory_generator.start_simulation(spatial_tessellation, mobility_diary, self.od_matrix)

            for time_slot, location in enumerate(synthetic_trajectory):
                f.write("%s,%s,%s\n" %(agent_id, location, time_slot))
                n_rows += 1

            # progress bar update
            percentage = int(float(count * 100) / self.n_agents)
            if percentage > old_p:
                hashes = '-' * int(round(percentage/5))
                spaces = ' ' * (bar_length - len(hashes))
                sys.stdout.write("\rExec: [{0}] {1}%".format(hashes + spaces, int(round(percentage))))
                sys.stdout.flush()
                old_p = percentage
            count += 1

        f.close()

        statinfo = os.stat(self.filename)
        sys.stdout.write("\n\nStats:\tfile lines %d\n\t\tfile size %s MB"
                         %(n_rows, float(statinfo.st_size) / (10 ** 6)))
        sys.stdout.flush()

if __name__ == '__main__':

    import argparse

    print "----------------------------------------------"
    print "                   DITRAS                 "
    print "        (DIary-based TRAjectory Simulator)   "
    print "----------------------------------------------"
    print "Authors: ", __author__
    print "Email:  ", __contact__
    print "----------------------------------------------\n"

    parser = argparse.ArgumentParser()
    parser.add_argument('n_agents', type=int, help='number of agents to simulate')
    parser.add_argument('length', type=int, help='length of the period in hours')
    parser.add_argument('spatial_tessellation', type=str, help='file where the spatial tessellation is stored')
    parser.add_argument('od_matrix', type=str, help='file where the od matrix is stored (if the file does not exists, '
                                                    'the od matrix will be computed and stored in a file with the '
                                                    'name specified in input)')
    parser.add_argument('diary_generator', type=str, help='the diary generator')
    parser.add_argument('filename', type=str, help='the name of the file where to store the trajectories')
    args = parser.parse_args()

    #max_locs = 1000
    #spatial_tessellation = dict(load_spatial_tessellation('location2info').items()[:max_locs])

    spatial_tessellation = load_spatial_tessellation(args.spatial_tessellation)  # de-commenta per matrice su tutte le locazioni
    if not path.isfile(args.od_matrix):
        compute_od_matrix(spatial_tessellation, filename=args.od_matrix)

    ditras = DITRAS(n_agents=args.n_agents, length=args.length, diary_generator=MD(filename=args.diary_generator),
                 trajectory_generator=dEPR(rho=0.6, gamma=0.21), filename=args.filename)

    ditras.start_simulation(spatial_tessellation, od_matrix=args.od_matrix)



