import numpy as np
from scipy.stats import poisson
from numpy import random
import seaborn as sb

class City:
    


    def __init__(self, city_name, population, transition_prob, days):
        """
        Constructor for City --- initilize city with variables:
        - String city_name: name of city / location modeling
        - double[] population: array of population size of each compartment correcponding to varaible compartments
        - String[] compartments: all the compartments that the disease progresses through
        - double [] policies: inter-compartmental transition probability 
        based on size of each compartment
            - policies[0]: transition prob from pre-symptomatic and asymptomatic to susceptible
            - policies[1]: transition probl from ill to susceptible
                - ill people are either hospitalized or isolated, therefore probability of 
                transitioning is much lower than those not showing symptom
        - dict() comparment_pop: dictionary that maps compartment to number of people in it
        - dict() daily_change: change in each comparment that occurs during a day
            - each distribution adds to the dict of daily change
        - days:  number of days we want to run the simulation for


        >>> san_diego = City("San Diego", [2000, 0, 0, 0, 0, 0, 0], [0.4, 0.1], 10)
        >>> san_diego.day_change_dist
        >>> print(san_diego.population_getter())
        [2000, 0, 0, 0, 0, 0, 0]
        >>> san_diego.transition()
        >>> san_diego.simulation()

        

        """
        self.city_name = city_name
        self.population = population

        self.compartments = ["Susceptible", "Exposed", "Pre-Symptomatic", "Asymptomatic",
                                "Ill", "Dead", "Recovered"]
        self.compartment_pop = dict(zip(self.compartments, self.population))
        self.days = range(0, 21)
        # 2d array as values of daily_change dict
        self.daily_change = [[0 for i in range(len(self.compartments))] for j in range(21)]
        self.day_change_dist = dict(zip(self.days, self.daily_change))
        

        self.transition_prob = transition_prob
       



    # getters for instance variables
    def city_getter(self):
        return self.city_name

    def population_getter(self):
        return self.population

    def compartment_pop_getter(self):
        return self.compartment_pop

    def transition_prob_getter(self):
        return self.transition_prob

    def day_change_dist(self):
        return self.day_change_dist



    # setter for instance variables
    def city_setter(self, new_city):
        self.city_name = new_city

    def population_setter(self, new_pop):
        self.population = new_pop
        self.compartment_pop = dict(zip(self.compartments, self.population))

    def transition_prob_setter(self, new_transition):
        self.transition_prob = new_transition




    # methods to alter transition_prob based on population for each compartment
    # for simulation purpose, these transitions will be combined into one method as one-day simulation
    
    # Susceptible ---> Exposed 
    def sus_expo_transition(self, rate, expo_prop):
        x = []
        for i in range(21):
            x.append(poisson.pmf(i, rate))
        for i in range(21):
            self.day_change_dist[i][1] += (x[i] * self.population[0] * expo_prop)
        

    # Exposed ---> Presymptomatic 
    def expo_presymp_transition(self, rate):
        x = []
        for i in range(21):
            x.append(poisson.pmf(i, rate))
        for i in range(21):
            self.day_change_dist[i][2] += (x[i] * self.population[1])
        


    # Presymptomatic ---> Asymptomatic Or Ill based on pre-symptomatic partition
    def presymp_asymp_illl_transition(self, rate, asmyp_prop):
        x = []
        for i in range(21):
            x.append(poisson.pmf(i, rate))
        for i in range(21):
            self.day_change_dist[i][3] += (x[i] * self.population[2] * asmyp_prop)
        y = []
        for i in range(21):
            y.append(poisson.pmf(i, rate))
        for i in range(21):
            self.day_change_dist[i][4] += (y[i] * self.population[2] * (1 - asmyp_prop))
        

    def ill_dead_recovered_transition(self, rate, death_prop):
        x = []
        for i in range(21):
            x.append(poisson.pmf(i, rate))
        for i in range(21):
            self.day_change_dist[i][5] += (x[i] * self.population[4] * death_prop)
        y = []
        for i in range(21):
            y.append(poisson.pmf(i, rate))
        for i in range(21):
            self.day_change_dist[i][6] += (x[i] * self.population[4] * (1 - death_prop))
        

    def transition(self):
        self.sus_expo_transition(2, 0.3)
        self.population[0] -= self.day_change_dist[0][1]
        self.population[1] += self.day_change_dist[0][1]
        self.expo_presymp_transition(3)
        self.population[1] -= self.day_change_dist[0][2]
        self.population[2] += self.day_change_dist[0][2]
        self.presymp_asymp_illl_transition(4, 0.8)
        self.population[2] -= self.day_change_dist[0][3]
        self.population[2] -= self.day_change_dist[0][4]
        self.population[3] += self.day_change_dist[0][3]
        self.population[4] += self.day_change_dist[0][4]
        self.ill_dead_recovered_transition(5, 0.2)
        self.population[4] -= self.day_change_dist[0][5]
        self.population[4] -= self.day_change_dist[0][6]
        self.population[5] += self.day_change_dist[0][5]
        self.population[6] += self.day_change_dist[0][6]


    def simulation(self):
        for i in range(self.days):
            self.transition()
            # print(self.day_change_dist)
            # print(self.population)
  







        

        

   
