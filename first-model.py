import numpy as np
from scipy.stats import poisson
from numpy import random
import matplotlib.pyplot as plt

class City:
    


    def __init__(self, city_name, population, transition_prob, days):
        """
        Constructor:
        - String city_name ---> city to be modeled
        - int[] population ---> population as an array of length 7: exisiting condition
        - double[] transition probability ---> transition probability of city, user can define based on own city condition
        - int num_days ---> numbers of dats to be simulated

        Attributes:
        - String[] compartments ---> all the stages of disease propogation, ordered in sequence of transmition
        - dict compartment_pop ---> dictionary with compartment as keys and population of each compartment as values
        - int[] days ---> from day 0 to day 20. used to distribute over Poisson distribution
        - dict day-change_dist ---> dictionary to record daily changes of population
            - looks something like this: {0: [s,e,p,a,i,d,r], 1: [s,e,p,a,i,d,r] ... 20: [s,e,p,a,i,d,r]}
        - int days_simulated ---> although desired simulation days is specified by num_days, there are conditions that will terminate the simulation. 
                                  this attribute will keep track of actual simulation days
        - boolean can_transmit ---> when False, terminate simulation
        - int[] susceptible ---> list to keep track of changes in susceptible compartment. used to plot population change



        >>> san_diego = City("San Diego", [1426000, 0, 0, 0, 0, 0, 0], [0.4, 0.1], 20)
        >>> print(san_diego.get_population())
        >>> san_diego.simulation()
        >>> print(san_diego.days_simulated)
        >>> san_diego.test_plot()


        

        """

        self.city_name = city_name
        self.population = population
        self.transition_prob = transition_prob
        self.num_days = days
        
        

        self.compartments = ["Susceptible", "Exposed", "Pre-Symptomatic", "Asymptomatic",
                                "Ill", "Dead", "Recovered"]
        self.compartment_pop = dict(zip(self.compartments, self.population))
        self.days = range(0, 21)
        self.daily_change = [[0 for i in range(len(self.compartments))] for j in range(21)]
        self.day_change_dist = dict(zip(self.days, self.daily_change))
        self.days_simulated = 0
        self.can_transmit = True
        self.susceptible = [self.population[0]]
        self.exposed = [self.population[1]]
        self.presymptomatic = [self.population[2]]
        self.asymptomatic = [self.population[3]]
        self.ill = [self.population[4]]
        self.dead = [self.population[5]]
        self.recovered = [self.population[6]]
       



    # getters for instance variables
    def get_city(self):
        return self.city_name

    def get_population(self):
        return self.population

    def get_compartment_pop(self):
        return self.compartment_pop

    def get_transition_prob(self):
        return self.transition_prob

    def get_day_change_dist(self):
        return self.day_change_dist

    

    # setter for instance variables
    def set_city(self, new_city):
        self.city_name = new_city

    def set_population(self, new_pop):
        self.population = new_pop
        self.compartment_pop = dict(zip(self.compartments, self.population))

    def set_transition_prob(self, new_transition):
        self.transition_prob = new_transition




    ## Transition between compartments ---> distribute case of transition into day0 to day20 using Poisson distribution    
    # Susceptible ---> Exposed 
    def sus_expo_transition(self, rate, expo_prop):
        # THE DAY OF self.transition() being called
        # expo_prop: proportion of susceptible population being exposed to virus
        x = []
        for i in range(21):
            x.append(poisson.pmf(i, rate))
            # if someone is susceptible, it might take from 0 to 20 days to become exposed, use Poisson distribution to model the probability of becoming exposed from each days
        for i in range(21):
            self.day_change_dist[i][1] += (self.population[0] * expo_prop * x[i])
            # of all the days from day0 to day20, the value correspond to change in exposed population in each values list gets updated
        

    # Exposed ---> Presymptomatic 
    def expo_presymp_transition(self, rate):
        # THE DAY OF self.transition() being called
        x = []
        for i in range(21):
            x.append(poisson.pmf(i, rate))
        for i in range(21):
            self.day_change_dist[i][2] += (x[i] * self.population[1])


    # Presymptomatic ---> Asymptomatic Or Ill based on pre-symptomatic partition
    def presymp_asymp_illl_transition(self, rate, asmyp_prop):
        # THE DAY OF self.transition() being called
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
        # THE DAY OF self.transition() being called
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
        self.sus_expo_transition(2, 0.03)
        self.expo_presymp_transition(3)
        self.presymp_asymp_illl_transition(4, 0.8)
        self.ill_dead_recovered_transition(5, 0.02)

        if(self.day_change_dist[0][1] <= self.population[0]):

            self.population[0] -= self.day_change_dist[0][1]
                    # susceptible population decrease
            self.population[1] += self.day_change_dist[0][1]
                    # exposed population increase
            self.population[1] -= self.day_change_dist[0][2]
                    # exposed population decrease
            self.population[2] += self.day_change_dist[0][2]
                    # presymptomatic population increase\
            self.population[2] -= self.day_change_dist[0][3]
                            # presymptomatic population decrease and become asymptomatic
            self.population[3] += self.day_change_dist[0][3]
                            # asymptomatic population increase
            self.population[3] -= self.day_change_dist[0][4]
                        # presymptomatic population decrease and become ill
            self.population[4] += self.day_change_dist[0][4]
                        # ill population increase
            self.population[4] -= self.day_change_dist[0][5]
                                    # ill population decrease to become dead
            self.population[5] += self.day_change_dist[0][5]
                                    # dead population increase
            self.population[5] -= self.day_change_dist[0][6]
                                    # ill population decrease to become recovered 
            self.population[6] += self.day_change_dist[0][6]
    
            self.susceptible.append(self.population[0])
            self.exposed.append(self.population[1])
            self.presymptomatic.append(self.population[2])
            self.asymptomatic.append(self.population[3])
            self.ill.append(self.population[4])
            self.dead.append(self.population[5])
            self.recovered.append(self.population[6])
            to_replace = list(self.day_change_dist.values())
            to_replace.pop(0)
            to_append = [0 for i in range(len(self.compartments))]
            to_replace.append(to_append)
            self.day_change_dist = dict(zip(self.days, to_replace))
                # print(self.day_change_dist)
        else:
            self.can_transmit = False






    def simulation(self):
        for i in range(self.num_days):
            self.transition()
            if(self.can_transmit == True): 
                self.days_simulated +=1;
            else:
                break;



    def test_plot(self):
        # print(self.exposed)
        # print(self.dead)
        x = range(self.days_simulated+1)
        plt.plot(x, self.susceptible, color = 'g')
        plt.plot(x, self.exposed, color = 'r')
        plt.plot(x, self.presymptomatic, color = 'm')
        plt.plot(x, self.asymptomatic, color = 'orange')
        plt.plot(x, self.ill, color = 'c')
        plt.plot(x, self.dead, color = 'y')
        plt.plot(x, self.recovered, color = 'b')
        plt.show()



            
  







        

        

   
