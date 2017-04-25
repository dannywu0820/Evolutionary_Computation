import pandas as pd
import sys
import matplotlib
import matplotlib.pyplot as plt
#import random
import numpy as np
from random import random, randint
from numpy import linalg
import itertools

def read_dataset(file_name):
    dataset = []
    fp = open(file_name, 'r')
    for line in fp.readlines():
        data = [float(line[0:7].strip()), float(line[7:16].strip())]
        dataset.append(data)

    return dataset

def plot_dataset(dataset):
    dataset = np.array(dataset)
    plt.scatter(dataset.T[0], dataset.T[1], color='black')
    plt.show()

def individual(length, min, max): #length is the length of gene
    return [ randint(min,max) for x in xrange(length) ]

def population(size, length, min, max):
    return [ individual(length, min, max) for x in xrange(size) ]

def fitness(individual, M):
    #print individual
    points_of_each_line = []
    mean_of_each_line = []
    direction_of_each_line = []
    MSD_of_each_line = [] #Mean Squared Distance of points
    for i in range(M):
        points_of_each_line.append(list())
        MSD_of_each_line.append(50000)

    #assign N points to M lines
    for i in range(len(individual)):
        points_of_each_line[individual[i]-1].append(dataset[i])

    for each in points_of_each_line:
        #print "-----[Line Info]-----"
        #print each
        if(len(each) > 1):
            #calculate mean of each set of points
            mean = np.mean(np.array(each), axis=0)
            mean_of_each_line.append(mean)
            #print mean

            #calculate covariance matrix of each set of points
            covariance_matrix = np.cov(np.array(each).T) #each row in transpose means a random variable
            #print covariance_matrix

            #find eigenvector of each covariance matrix, the largest eigenvector is the direction
            eig_val, eig_vec = linalg.eig(covariance_matrix) #can test with np.array([[1,2],[4,3]])
            if(eig_val[0] > eig_val[1]):
                direction_of_each_line.append(eig_vec[0])
            else:
                direction_of_each_line.append(eig_vec[1])
            #print eig_val
            #print eig_vec
        else:
            mean_of_each_line.append(np.array([0,0]))
            direction_of_each_line.append(np.array([0,0]))

    #calculate MSD of points to their assigned lines
    for i in range(M):
        num_of_points = len(points_of_each_line[i])
        if(num_of_points > 1):
            sqr_dis = 0
            for point in points_of_each_line[i]:
                dis = distance_to_line(np.array(point), mean_of_each_line[i], direction_of_each_line[i])
                sqr_dis+=(dis*dis)
            MSD = sqr_dis/num_of_points
            MSD_of_each_line[i] = MSD

    #print MSD_of_each_line
    fitness_value = sum(MSD_of_each_line)/M #1/sum(MSD_of_each_line)
    #print "fitness value = "+str(fitness_value)
    return fitness_value

def distance_to_line(point, mean, direction):
    perpendicular_direction = np.array([direction[1], -direction[0]]) #the vector which is perpendicular to direction vector
    PQ = point-mean
    distance = linalg.norm(np.dot(PQ, perpendicular_direction))/linalg.norm(perpendicular_direction)
    return distance

def crossover(mate1, mate2, crossover_type, crossover_rate):
    child = []
    if crossover_type == '1-point':
        crossover_point = randint(0, len(mate1)-1)
        child = mate1[:crossover_point] + mate2[crossover_point:]
    elif crossover_type == '2-point':
        crossover_point = randint(0, len(mate1)-1)
        crossover_point2 = randint(0, len(mate2)-1)
        if(crossover_point < crossover_point2):
            child = mate1[:crossover_point] + mate2[crossover_point:crossover_point2] + mate1[crossover_point2:]
        else:
            child = mate1[:crossover_point2] + mate2[crossover_point2:crossover_point] + mate1[crossover_point:]
    elif crossover_type == 'uniform':
        for i in range(0, len(mate1)):
            if crossover_rate > random():
                child.append(mate1[i])
            else:
                child.append(mate2[i])
    else:
        pass

    return child

def mutate(individual, mutate_type, mutate_rate):
    individual_len = len(individual)
    if mutate_type == 'swap':
        if mutate_rate > random():
            mutate_point = randint(0, individual_len-1)
            mutate_point2 = randint(0, individual_len-1)
            individual[mutate_point], individual[mutate_point2] = individual[mutate_point2], individual[mutate_point]
    elif mutate_type == 'inversion':
        if mutate_rate > random():
            mutate_point = randint(0, individual_len-1)
            mutate_point2 = randint(mutate_point, individual_len-1)
            reversed_part = list(reversed(individual[mutate_point:mutate_point2]))
            #print individual[:mutate_point]
            #print reversed_part
            #print individual[mutate_point2:]
            individual = individual[:mutate_point] + reversed_part + individual[mutate_point2:]
    else:
        if mutate_rate > random():
            mutate_point = randint(0, individual_len-1)
            individual[mutate_point] = randint(1, 4)

    return individual

evolve_history = []
def evolve(pop, groups, control_flow_type='steady-state', retain_rate=0.01, crossover_rate=0.5, mutate_rate=1, random_select=0.01):
    grade = [ (fitness(ind, groups), ind) for ind in pop]
    #sorted in ascending order, the lower a fitness value is, the better the individual is
    sorted_pop = [ x[1] for x in sorted(grade) ]
    print fitness(sorted_pop[0], groups)
    evolve_history.append(fitness(sorted_pop[0], groups))

    retained_len = int(len(sorted_pop)*retain_rate)
    parents = sorted_pop[:retained_len] #form mating pools

    #add other individuals to promote genetic diversity
    for individual in sorted_pop[retained_len:]:
        if random_select > random():
            parents.append(individual)

    #crossover parents to create children
    individual_len = len(pop[0])
    parents_len = len(parents)
    desired_len = len(pop) - parents_len
    children = []
    if control_flow_type != 'generational':
        #keep best parents in next generation
        children.extend(parents)
    while len(children) < len(pop):
        mate1 = randint(0, parents_len-1)
        mate2 = randint(0, parents_len-1)
        mate1 = parents[mate1]
        mate2 = parents[mate2]
        crossover_point = randint(0, individual_len-1)
        child = crossover(mate1, mate2, '1-point', crossover_rate)
        #child = mate1[:crossover_point] + mate2[crossover_point:]
        #print str(child.count(1)) + ":" + str(child.count(2)) + ":" + str(child.count(3)) + ":" + str(child.count(4))
        children.append(child)

    #mutate some individuals
    for individual in children:
        individual = mutate(individual, 'swap', mutate_rate)

    return children #parents

def genetic_algorithm(pop_size=50, ind_len=100, gene_min=1, gene_max=4, gen_num=200):
    pop = population(pop_size, ind_len, gene_min, gene_max)
    for i in range(0, gen_num):
        pop = evolve(pop, groups=gene_max)

    final_grade = [ (fitness(ind, gene_max), ind) for ind in pop]
    best_ind = sorted(final_grade)[0][1]
    print best_ind

    plot_result(best_ind)

def plot_result(individual):
    colors = ['red', 'blue', 'green', 'black', 'magenta']
    for gene,xy in zip(individual, dataset):
        plt.scatter(xy[0], xy[1], color=colors[gene-1])
    plt.show()

    evolve_gen = list(range(1, len(evolve_history)+1))
    plt.plot(evolve_gen, evolve_history, color='b')
    plt.show()

def get_permutations(range=0, repeat=2):
    x = [1,2,3]
    all_permutations = [p for p in itertools.product(x, repeat=repeat)]
    print all_permutations

if __name__ == '__main__':
    dataset = read_dataset('./Dataset/lineN100M4.txt')
    #plot_dataset(dataset)
    genetic_algorithm(pop_size=100, ind_len=100, gene_min=1, gene_max=4, gen_num=200)
    #get_permutations()
