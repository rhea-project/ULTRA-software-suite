#-----------------------------------------------------------------------------#
# Import modules
#-----------------------------------------------------------------------------#
# standard library imports
import numpy as np

# related third party imports
from tqdm import tqdm
#import neurolab as nl

# local application/library specific imports

#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#class NeuralNetwork:
#    def __init__(self):
#        self.nodeLayer = 10
#        self.input      = None
#        self.y          = None
#        
#    def reShape(self):
#        self.weights1   = 2*np.random.rand(self.input.shape[1],self.nodeLayer)-1
#        self.weights2   = 2*np.random.rand(self.nodeLayer,1)-1
#        self.output     = np.zeros(self.y.shape)
#        
#    def sigmoid(self, x):
#        return 1 / (1 + np.exp(-x))
#    
#    def sigmoid_derivative(self, x):
#        return x * (1 - x)
#        
#    def feedforward(self):
#        self.layer1 = self.sigmoid(np.dot(self.input, self.weights1))
#        self.output = self.sigmoid(np.dot(self.layer1, self.weights2))
#
#    def backprop(self):
#        d_weights2 = np.dot(self.layer1.T, 
#                            (2*(self.y - self.output) * self.sigmoid_derivative(self.output)))
#        d_weights1 = np.dot(self.input.T, 
#                            (np.dot(2*(self.y - self.output) * self.sigmoid_derivative(self.output), self.weights2.T) * self.sigmoid_derivative(self.layer1)))
#
#        # update the weights with the derivative (slope) of the loss function
#        self.weights1 += d_weights1
#        self.weights2 += d_weights2
#        
#    def train(self, nb_iter):
#        for iter in range(nb_iter):
#            self.feedforward()
#            self.backprop()
##            print(self.y - self.output)
    
        
        
class DifferentialEvolution(object):
    def __init__(self, fobj, bounds, mut=0.8, crossp=0.7, popsize=8, its=1, mod=None):
        
        self.x_data_norm = []
        self.f_data_norm = [] 
        
        self.x_data_denorm = []
        self.f_data_denorm = [] 
        
        self.f_data_max = 0
        
        self.bounds = bounds
        self.mut = mut
        self.crossp = crossp
        self.popsize = popsize
        self.its = its
               
#        self.metaModell = nl.net.newff([ [0,1] for bt in bounds ], [10, 1])
        
        self.fobj = fobj
    
    def meta(self):
#        inp = np.asarray(self.x_data_list[0]) / 1e+08
#        temp = np.asarray(self.f_data_list[0])
#        tar = np.reshape(temp, (len(temp), 1)) / 1e+11
#        error = self.metaModell.train(inp, tar, epochs=500, show=100, goal=0.02)    
        inp = np.asarray(self.x_data_norm)
        tar = np.asarray(self.f_data_norm)
#        tar = np.reshape(temp, (len(temp), 1))
        error = self.metaModell.train(inp, tar, epochs=500, goal=0.001)
    
#    def normalizeData(self):
##        self.f_data_norm = list(self.f_data_denorm / max(self.f_data_denorm)
#        self.f_data_norm = []
#        self.f_data_max = np.asarray(max(self.f_data_denorm))
#        for f in range(len(self.f_data_denorm)):
#            temp = np.asarray(self.f_data_denorm[f]) / self.f_data_max
#            self.f_data_norm.append(list(temp))
            
    def normX(self, x):
        min_b, max_b = np.asarray(self.bounds).T
        diff = np.fabs(min_b - max_b)
        x_norm = (x - min_b) / diff
        
        return x_norm
    
    def denormData(self, x):
        x_norm = self.min_f + x * self.diff
        f_norm = np.exp(x_norm)
        
        return f_norm
        
    
    def normalizeData(self):
        self.min_f = 0
        self.max_f = 0
        self.diff = 0
        self.f_data_norm = []
        x = np.asarray(self.f_data_denorm)
        x_denorm_log = np.log(x)
        
        min_f = np.min(x_denorm_log)
        max_f = np.max(x_denorm_log)
                
        diff = max_f - min_f
        
        x_norm_log = ( x_denorm_log - min_f ) / diff
        
        self.min_f = min_f
        self.max_f = max_f
        self.diff = diff
        self.f_data_norm = list(x_norm_log)
    
    def denormX(self, x):
        min_b, max_b = np.asarray(self.bounds).T
        diff = np.fabs(min_b - max_b)
        x_denorm = min_b + x * diff
        x_log = np.log(x_denorm)
        
        min_f = np.min(x_log)
        max_f = np.max(x_log)
        
        diff = max_f - min_f
        
        x_denorm_log = min_f + x * diff
        
        return x_denorm
    
    def initialPopulation(self):
        # dimension of the problem: number of design variables
        dimensions = len(self.bounds)
        
        # random initial population 
        # TODO: systematic creation of initial population (DOE)
        pop = np.random.rand(self.popsize, dimensions)
        
        # denormalizing x obtained from initial population
        min_b, max_b = np.asarray(self.bounds).T
        diff = np.fabs(min_b - max_b)
        pop_denorm = min_b + pop * diff
        
        
        # evaluate fobj with initial population
        fitness = np.asarray([self.fobj(ind) for ind in pop_denorm])
        
        for p in range(len(pop)):
            self.x_data_norm.append(list(pop[p]))
            self.x_data_denorm.append(list(pop_denorm[p]))
            self.f_data_denorm.append(list(fitness[p]))
        
        # find index of fittest initial population
        best_idx = np.argmin(fitness)
        
        # fitness of fittest popolation
        best = pop_denorm[best_idx]
        
        self.min_b = min_b
        self.max_b = max_b
        
        self.diff = diff
        self.pop = pop
        
        self.dimensions = dimensions
        
        self.fitness = fitness
        self.best_idx = best_idx
        self.best = best
        
        
    def de(self):
            """ Initialization that builds the port.
            
                Assumptions:
                None
            
                Source:
                N/A
            
                Inputs:
                fobj
                bounds
            
                Outputs:
                N/A
                
                Properties Used:
                N/A
            """          
                        
            # dimension of the problem: number of design variables
            dimensions = len(self.bounds)
            
            # random initial population 
            # TODO: systematic creation of initial population (DOE)
            pop = np.random.rand(self.popsize, dimensions)
            
            # denormalizing x obtained from initial population
            min_b, max_b = np.asarray(self.bounds).T
            diff = np.fabs(min_b - max_b)
            pop_denorm = min_b + pop * diff            
            
            # evaluate fobj with initial population
            fitness = np.asarray([self.fobj(ind) for ind in pop_denorm])
            
#            for p in range(len(pop)):
#                self.x_data_norm.append(list(pop[p]))
#                self.x_data_denorm.append(list(pop_denorm[p]))
#                self.f_data_denorm.append(list(fitness[p]))
            
            # find index of fittest initial population
            best_idx = np.argmin(fitness)
            
            # fitness of fittest popolation
            best = pop_denorm[best_idx]
            
#            self.min_b = min_b
#            self.max_b = max_b
#            
#            self.diff = diff
#            self.pop = pop
#            
#            self.dimensions = dimensions
#            
#            self.fitness = fitness
#            self.best_idx = best_idx
#            self.best = best
        
            # iterate over number of iterations
            for i in tqdm(range(self.its)):
                
#                self.normalizeData()
#                self.meta()
                
                # iterate over initial population
                for j in range(self.popsize):
                    # generate list with indices of population vectors excluding 
                    # target vector
                    idxs = [idx for idx in range(self.popsize) if idx != j]
                    
                    # randomly choose 3 indices
                    selected = np.random.choice(idxs, 3, replace=False)
                    
                    # get population corresponing to randomly chosen indices
                    a, b, c = pop[selected]
                    
                    # create mutant vector by computing difference of b and c 
                    # and add this difference weighted by mutation factor mut to a
                    # mut usally in [0.5, 2.0] --> the higher the larger search radius
                    # but slows down convergence
                    mutant = a + self.mut * (b - c)
                    
                    # clippling values of mutant to 0 and 1 
                    mutant = np.clip(mutant, 0, 1)
                    
                    # recombination
                    # mixing information of current vector and mutant
                    # decide which indices selected from mutant
                    cross_points = np.random.rand(dimensions) < self.crossp
                    if not np.any(cross_points):
                        cross_points[np.random.randint(0, dimensions)] = True
                    
                    # recombine mutant and current vector to trial vector
                    trial = np.where(cross_points, mutant, pop[j])
                    
                    # denormalize trial vector
                    trial_denorm = min_b + trial * diff
        
#                    self.x_data_norm.append(list(trial))
#                    self.x_data_denorm.append(list(trial_denorm))
                    
                    # evaluate trial vector
                    f = self.fobj(trial_denorm)
        
#                    self.f_data_denorm.append(list(f))
                    
                    # compare fitness of trial to current vector and refresh fittest
                    # vector
                    if f < fitness[j]:
                        fitness[j] = f
                        pop[j] = trial
                        if f < fitness[best_idx]:
                            best_idx = j
                            best = trial_denorm 
                                                
                yield best, fitness[best_idx]

# 
#def fobj(x):
#    return x[0]**2 + x[1]**3 + x[2]
#
#bounds = [(0,10), (0,10), (0,10)]
#DE = DifferentialEvolution(fobj, bounds)
#DE.initialPopulation()
#
#
#    
#NN = NeuralNetwork(DE.x_data, DE.f_data)
##NN = NeuralNetwork(x, y)
#NN.training(3)
## testing
#NN.input = np.array([5.00802414, 1.08458117, 7.44184035])
#NN.feedforward()
#print(NN.output)

#def de(self, fobj, bounds, mut=0.8, crossp=0.7, popsize=8, its=10, mod=None):
#    """ Initialization that builds the port.
#    
#        Assumptions:
#        None
#    
#        Source:
#        N/A
#    
#        Inputs:
#        fobj
#        bounds
#    
#        Outputs:
#        N/A
#        
#        Properties Used:
#        N/A
#    """          
#    
#    # dimension of the problem: number of design variables
#    dimensions = len(bounds)
#    
#    # random initial population 
#    # TODO: systematic creation of initial population (DOE)
#    pop = np.random.rand(popsize, dimensions)
#    
#    # denormalizing x obtained from initial population
#    min_b, max_b = np.asarray(bounds).T
#    diff = np.fabs(min_b - max_b)
#    pop_denorm = min_b + pop * diff
#    
#    # evaluate fobj with initial population
#    fitness = np.asarray([fobj(ind) for ind in pop_denorm])
#    
#    # find index of fittest initial population
#    best_idx = np.argmin(fitness)
#    
#    # fitness of fittest popolation
#    best = pop_denorm[best_idx]
#    
#    # iterate over number of iterations
#    for i in tqdm(range(its)):
#        
#        # iterate over initial population
#        for j in range(popsize):
#            # generate list with indices of population vectors excluding 
#            # target vector
#            idxs = [idx for idx in range(popsize) if idx != j]
#            
#            # randomly choose 3 indices
#            selected = np.random.choice(idxs, 3, replace=False)
#            
#            # get population corresponing to randomly chosen indices
#            a, b, c = pop[selected]
#            
#            # create mutant vector by computing difference of b and c 
#            # and add this difference weighted by mutation factor mut to a
#            # mut usally in [0.5, 2.0] --> the higher the larger search radius
#            # but slows down convergence
#            mutant = a + mut * (b - c)
#            
#            # clippling values of mutant to 0 and 1 
#            mutant = np.clip(mutant, 0, 1)
#            
#            # recombination
#            # mixing information of current vector and mutant
#            # decide which indices selected from mutant
#            cross_points = np.random.rand(dimensions) < crossp
#            if not np.any(cross_points):
#                cross_points[np.random.randint(0, dimensions)] = True
#            
#            # recombine mutant and current vector to trial vector
#            trial = np.where(cross_points, mutant, pop[j])
#            
#            # denormalize trial vector
#            trial_denorm = min_b + trial * diff
#            
#            # evaluate trial vector
#            f = fobj(trial_denorm)
#            
#            # compare fitness of trial to current vector and refresh fittest
#            # vector
#            if f < fitness[j]:
#                fitness[j] = f
#                pop[j] = trial
#                if f < fitness[best_idx]:
#                    best_idx = j
#                    best = trial_denorm 
#                                        
#        yield best, fitness[best_idx]
