import pandas as pd
import sys
import matplotlib
import matplotlib.pyplot as plt
#import random
import numpy as np
from random import random, randint
from numpy import linalg

#Step1.Read Data
dataset = []
file_location = './Dataset/lineN200M3.txt'
fp = open(file_location, 'r')
for line in fp.readlines():
    #print line[0:7].strip() + ',' + line[7:16].strip()
    data = [float(line[0:7].strip()), float(line[7:16].strip())]
    dataset.append(data)

#dataset = np.array(dataset)
#print population
#plt.scatter(dataset.T[0], dataset.T[1], color='b')
#plt.show()

def individual(length, min, max): #length is the length of gene
    return [ randint(min,max) for x in xrange(length) ]
#print individual(5, 1, 3)

def population(size, length, min, max):
    return [ individual(length, min, max) for x in xrange(size) ]

def fitness(individual, M):
    #print individual
    points_in_each_line = []
    mean_of_each_line = []
    direction_of_each_line = []
    MSD_of_each_line = []
    for i in range(M):
        points_in_each_line.append(list())
        MSD_of_each_line.append(50000)
    #divide into M sets
    for i in range(len(individual)):
        points_in_each_line[individual[i]-1].append(dataset[i])

    for each in points_in_each_line:
        #print "-----[Line]-----"
        #print each
        if(len(each) != 0 and len(each) != 1):
            #calculate mean of each set of points
            mean_of_each_line.append(np.mean(np.array(each), axis=0))
            #print mean_of_each_line

            #calculate covariance matrix of each set
            #print np.array(each).T #each row in transpose means a random variable
            covariance_matrix = np.cov(np.array(each).T)
            #print covariance_matrix

            #find eigenvector of each covariance matrix
            #eig_val, eig_vec = linalg.eig(np.array([[1,2],[4,3]]))
            eig_val, eig_vec = linalg.eig(covariance_matrix)
            #print eig_val
            #print eig_vec
            if(eig_val[0] > eig_val[1]):
                direction_of_each_line.append(eig_vec[0])
            else:
                direction_of_each_line.append(eig_vec[1])
        else:
            mean_of_each_line.append(np.array([0,0]))
            direction_of_each_line.append(np.array([0,0]))


    #calculate MSD of each line
    for i in range(M):
        num_of_points = len(points_in_each_line[i])
        if(num_of_points > 1):
            sqr_dis = 0
            for point in points_in_each_line[i]:
                dis = distance(point, mean_of_each_line[i], direction_of_each_line[i])
                sqr_dis+=(dis*dis)
            MSD = sqr_dis/num_of_points
            MSD_of_each_line[i] = MSD

    #for i in range(M):
        #print MSD_of_each_line[i]
    fitness_value = sum(MSD_of_each_line) #1/sum(MSD_of_each_line)
    #print "fitness value = "+str(fitness_value)
    return fitness_value

def distance(point, mean, direction):
    point_out = np.array(point)
    pointB = mean
    pointC = mean + direction*10
    #distance = linalg.norm(np.cross(pointB-pointC, pointC-point_out))/linalg.norm(pointB-pointC)
    #print distance

    perpen_dir = np.array([direction[1], -direction[0]]) #the vector which is perpendicular to direction
    PQ = point_out-pointB
    distance = linalg.norm(np.dot(PQ, perpen_dir))/linalg.norm(perpen_dir)
    return distance

#fitness(individual(100, 1, 3), 3)

evolve_history = []
def evolve(pop, groups, retain=0.2, crossover=0.05, mutate=0.01):
    grade = [ (fitness(ind, groups), ind) for ind in pop]
    #sorted in ascending order, the lower fitness value is, the better individual is
    sorted_pop = [ x[1] for x in sorted(grade) ]
    print fitness(sorted_pop[0], groups)
    evolve_history.append(fitness(sorted_pop[0], groups))
    #print "[Current Generation]"
    #print np.array(sorted_pop)

    retained_len = int(len(sorted_pop)*retain)
    parents = sorted_pop[:retained_len] #form mating pools
    #print parents

    #add other individuals to promote genetic diversity
    #crossover parents to create children
    individual_len = len(pop[0])
    parents_len = len(parents)
    desired_len = len(pop) - parents_len
    children = []
    while len(children) < desired_len:
        mate1 = randint(0, parents_len-1)
        mate2 = randint(0, parents_len-1)
        mate1 = parents[mate1]
        mate2 = parents[mate2]
        crossover_point = randint(0, individual_len-1)
        child = mate1[:crossover_point] + mate2[crossover_point:]
        children.append(child)
    parents.extend(children)
    #print "[Next Generation]"
    #print np.array(parents)
    #mutate some individuals
    for ind in parents:
        #for i in range(individual_len):
        if random() > mutate:
            #print "mutate"
            mutate_pos = randint(0, individual_len-1)
            ind[mutate_pos] = randint(1, groups)

    return parents

def plot_result(individual):
    colors = ['red', 'blue', 'green', 'black', 'magenta']
    for gene,xy in zip(individual, dataset):
        plt.scatter(xy[0], xy[1], color=colors[gene-1])
    plt.show()

    plt.plot(evolve_gen, evolve_history, color='b')
    plt.show()

pop_size = 100
ind_len = 200
ind_min = 1
ind_max = 3
pop = population(pop_size, ind_len, ind_min, ind_max)

gen_num = 200
evolve_gen = list(range(0, gen_num))
for i in range(0, gen_num):
    pop = evolve(pop, groups=ind_max)

final_grade = [ (fitness(ind, ind_max), ind) for ind in pop]
best_ind = sorted(final_grade)[0][1]
print best_ind

plot_result(best_ind)
