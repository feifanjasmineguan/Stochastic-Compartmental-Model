import numpy as np
from scipy.stats import poisson
from numpy import random
import seaborn as sb
import matplotlib.pyplot as plt

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


        >>> san_diego = City("San Diego", [1426000, 0, 0, 0, 0, 0, 0], [0.4, 0.1], 20)
        >>> print(san_diego.population_getter())
        >>> san_diego.simulation()
        >>> print(san_diego.days_simulated)


        >>> san_diego.test_plot()


        

        """
        self.num_days = days
        self.city_name = city_name
        self.population = population

        self.compartments = ["Susceptible", "Exposed", "Pre-Symptomatic", "Asymptomatic",
                                "Ill", "Dead", "Recovered"]
        self.compartment_pop = dict(zip(self.compartments, self.population))
        self.days = range(0, 21)
        # 2d array as values of daily_change dict
        self.daily_change = [[0 for i in range(len(self.compartments))] for j in range(21)]
        self.day_change_dist = dict(zip(self.days, self.daily_change))
        self.days_simulated = 0

        self.transition_prob = transition_prob
        self.can_transmit = True
        self.susceptible = [self.population[0]]
        self.exposed = [self.population[1]]
        self.presymptomatic = [self.population[2]]
        self.asymptomatic = [self.population[3]]
        self.ill = [self.population[4]]
        self.dead = [self.population[5]]
        self.recovered = [self.population[6]]
       



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
        # THE DAY OF self.transition() being called
        x = []
        for i in range(21):
            x.append(poisson.pmf(i, rate))
        for i in range(21):
            self.day_change_dist[i][1] += (x[i] * self.population[0] * expo_prop)
        

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
        if(self.day_change_dist[0][1] <= self.population[0]):

            self.population[0] -= self.day_change_dist[0][1]
                    # susceptible population decrease
            self.population[1] += self.day_change_dist[0][1]
                    # exposed population increase
    
            self.expo_presymp_transition(3)
            if(self.day_change_dist[0][2] <= self.population[1]):
                self.population[1] -= self.day_change_dist[0][2]
                            # exposed population decrease
                self.population[2] += self.day_change_dist[0][2]
                            # presymptomatic population increase
            else:
                self.can_transmit = False
                return
            
            self.presymp_asymp_illl_transition(4, 0.8)
            if(self.day_change_dist[0][3] <= self.population[2] and self.day_change_dist[0][4] <= self.population[3]):
                self.population[2] -= self.day_change_dist[0][3]
                        # presymptomatic population decrease and become asymptomatic
                self.population[3] += self.day_change_dist[0][3]
                        # asymptomatic population increase
                self.population[3] -= self.day_change_dist[0][4]
                    # presymptomatic population decrease and become ill
                self.population[4] += self.day_change_dist[0][4]
                    # ill population increase
            else:
                self.can_transmit = False
                return

            self.ill_dead_recovered_transition(5, 0.02)
            if(self.day_change_dist[0][5] <= self.population[4] and self.day_change_dist[0][6] <= self.population[5]):
                self.population[4] -= self.day_change_dist[0][5]
                        # ill population decrease to become dead
                self.population[5] += self.day_change_dist[0][5]
                        # dead population increase
                self.population[5] -= self.day_change_dist[0][6]
                        # ill population decrease to become recovered 
                self.population[6] += self.day_change_dist[0][6]
                        # recovered population increase
            else:
                self.can_transmit = False
                return 

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



            
  







        

        

   
