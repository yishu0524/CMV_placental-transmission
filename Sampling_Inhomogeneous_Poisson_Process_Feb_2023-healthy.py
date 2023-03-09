#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import math
from matplotlib import pyplot as plt
import random
import pandas as pd

import seaborn as sns




def sample_from_exponential_dist(a):
    # random sampling from exponential distribution with mean a
    # we will use inverse sampling 
    u = random.random()
    x = (1/a) * math.log(1/(1-u))
    return x




def inhomogeneousPoisson(Lambda, T):
    '''
    First, we aim to generate an array arrival_times[] representing the arrival time from an homogeneous Poisson Process with 
    fixed rate lambda_max up to time T; meanwhile, we also generate samples from uniform distribution on [0,1] which will decide
    whether these arrival times from homogeneous Poisson Process will be accepted as sample of the inhomogeneous Poisson Process
    '''
    lambda_max = max(Lambda[i][1] for i in range(len(Lambda)))
    
    arrival_times = [0]
    unif_samples = []
    while arrival_times[-1] <= T:
        arrival_times.append(arrival_times[-1]+sample_from_exponential_dist(lambda_max))
        unif_samples.append(random.random())
        
    # take out the first dummy time 0 and the last time which exceeds maximum time T     
    arrival_times = arrival_times[1:-1]
    
    # check if we accept each arrival time in arrival_times[]
    accepted_arrival_times = []
    for i in range(len(arrival_times)):
        # find the rate lambda[t] at the proposed arrival time
        for j in range(len(Lambda)):
            if Lambda[j][0] > arrival_times[i]:
                break
        acceptance_probability = Lambda[j][1]/lambda_max
        if acceptance_probability > unif_samples[i]:
            accepted_arrival_times.append(arrival_times[i])
    
    return accepted_arrival_times


# In[2]:


class CMVsystem:

    def __init__(self, lam, mu, beta, alpha, delta,S0,p,c, S, L, I, V, ratio,virus_entry,threshold):
        self.lam = lam
        self.mu = mu
        self.beta = beta
        self.alpha = alpha
        self.delta=delta
        self.p = p
        self.c = c
        self.S0= S0
        self.t = 0.
        self.S = S
        self.L = L
        self.I = I
        self.V = V
        self.virus_entry=virus_entry
        self.ratio=ratio
        self.trajectory = np.array([[self.S, self.L, self.I,self.V]],dtype=float)
        self.times = None
        self.threshold = threshold
            
#Define an initializer (__init__(self,beta,gamma,S,I,R)) which accepts model parameters beta and gamma, 
#and initial numbers of hosts in each of the S,I and R compartments.
    def reset(self, S, I, L, V,t=0.):
        self.t = t
        self.S = S
        self.L = L
        self.I = I
        self.V = V
        self.trajectory = np.array([[self.S, self.L, self.I, self.V]],dtype=float)


# In[3]:


class StochasticCMVsystem (CMVsystem):

    """Define a specialized subclass of the general SIRsystem for modeling SIR dynamics as a stochastic, continuous
    time process, using the Gillespie method of continuous time Monte Carlo"""

    def step(self):
        """Implement one step of Gillespie's Direct Method based on current reaction rates: identify a reaction to fire
        next, as well as a time at which it will fire, and execute that reaction (similar to as was described in the
        StochasticCells exercise)."""
        
        #if no virus inside the system, wait till the first virus enters the system
        if self.V + self.I + self.L == 0:
            if self.virus_entry:
                self.t = virus_entry.pop()
                self.V = 1
            else:
                return None, self.t
        
            
        
        
        s = np.random.uniform(0,1)
        transition = None
        #S_birth = self.lam
        #S_death = self.mu*self.S
        inf_rate = self.beta*self.V*self.S
        L_to_I = self.alpha*self.L
        #L_death = self.mu*self.L
        I_death = self.delta*self.I
        V_production = self.p*self.I
        V_death = self.c*self.V
        
        
        #S is kept constant
        #CASE1:with latent cells
        
        if s > self.ratio:
            rates=[inf_rate,L_to_I,I_death,V_production,V_death]
            total_rate=np.sum(rates)
            if total_rate == 0.:
                dt=0
            else:
                dt = np.random.exponential(1./total_rate, 1)[0]
            
            if self.virus_entry and self.t+dt > self.virus_entry[-1]:
                #print(self.t , self.virus_entry[-1])
                #recation does not happen, virus entry happens(case 1)
                self.V += 1
                self.t = self.virus_entry.pop()
                transition = 0
                
            else:
                #reaction happens
                if total_rate == 0.:
                    return transition, self.t
                if self.V+self.I+self.L == 0:
                    return transition, self.t
                if self.V+self.I+self.L >= self.threshold:
                    return transition, self.t
                ranno = np.random.uniform(0,1)
                if ranno < np.sum(rates[0:1])/total_rate:
                    self.L += 1   #Latent cell production
                    transition = 1
                elif ranno < np.sum(rates[0:2])/total_rate:
                    self.I+=1     #Latent cell become infected
                    self.L-=1
                    transition = 2
                elif ranno < np.sum(rates[0:4])/total_rate:
                    self.I-=1      #Infected Cell death
                    transition = 3
                elif ranno < np.sum(rates[0:5])/total_rate:
                    self.V+=1      #virus creation
                    transition = 4
                elif ranno < np.sum(rates[0:6])/total_rate:
                    self.V-=1     #virus death
                    transition = 5
                    
                self.t += dt
                

        if s <= self.ratio:
            
            rates=[inf_rate,I_death,V_production,V_death]
            total_rate=np.sum(rates)
            if total_rate == 0.:
                dt=0
            else:
                dt = np.random.exponential(1./total_rate, 1)[0]
                

            if self.virus_entry and self.t+dt > self.virus_entry[-1]:
                #recation does not happen, virus entry happens
                self.V += 1
                self.t = self.virus_entry.pop()
                transition = 0
                
            else:
                #reaction happens

                if total_rate == 0.:
                    return transition, self.t
                if self.V + self.I==0:
                    return transition, self.t
                if self.V + self.I>=self.threshold:
                    return transition, self.t

                ranno = np.random.uniform(0,1)        

                if ranno < np.sum(rates[0:1])/total_rate:
                    self.I += 1
                    transition = 1

                elif ranno < np.sum(rates[0:2])/total_rate:
                    self.I -= 1
                    transition = 2
                elif ranno < np.sum(rates[0:3])/total_rate:
                    self.V+=1
                    transition = 3
                elif ranno < np.sum(rates[0:4])/total_rate:
                    self.V-=1
                    transition = 4
                self.t += dt
                
                
                
        #dt = np.random.exponential(1./total_rate, 1)[0]
        #why 1/the hazard function? different from the textbook
        #print(transition, self.t)
        return transition, self.t

            
            
    def run(self, T = None, make_traj=True):
        """Run the Gillespie algorithm for stochastic simulation from time 0 to at least time T, starting with the initial
        values stored in the S,I,R state variables; story the result in self.trajectory if make_traj argument is
        set to True"""

        if T is None:
            T = sys.maxsize
        self.times = [0.]

        transition = 1
        while self.t < T:
            transition, t = self.step()
            if transition == None:
                return self.t
            if make_traj:
                self.trajectory = np.concatenate((self.trajectory, [[self.S,self.L,self.I,self.V]]), axis=0)
            self.times.append(self.t)
        return self.t


# In[4]:


count=0;


# In[5]:


infected_time=[];
infected_status=[];


# In[8]:


healthy_pop=pd.read_csv('healthy populationflow.csv')


# In[32]:


infected_time=[];
infected_status=[];
for index in range(1000):
    
    Lambda= list(zip(healthy_pop['days'],healthy_pop['pop_flow']))
    arrival_times = inhomogeneousPoisson(Lambda, 166)
    virus_entry=arrival_times[::-1]
    threshold = 800
    mu=1/4.5; beta=3e-12; alpha=1; delta=0.77; S0=4e8;lam=mu*S0 ;p=1600;c=2; S=S0; L=0; I=0; V=0;ratio=0.5;
    M1 = StochasticCMVsystem(lam, mu, beta, alpha, delta, S0, p, c, S, L, I, V,ratio,virus_entry,threshold);
    
    M1.run(166)
    infected_time.append(M1.times[:-1][-1])
    if M1.trajectory[-1,3]>=threshold:
        count+=1
        infected_status.append(1)
    else:
        infected_status.append(0)
    if index%50 == 0:
        print(index)

    


# In[33]:


duration=[]

for i,x in enumerate(infected_status):
    if x >0:
        duration.append(infected_time[i])
    else:
        duration.append(166)

d3 = {'duration': duration,'event':infected_status}
df4 = pd.DataFrame(data=d3)

#df4.to_csv('health_events_500_shift_00.csv')


# In[ ]:


# Python code to create the above Kaplan Meier curve
from lifelines import KaplanMeierFitter

## Example Data 
durations = df4['duration']
event_observed = df4['event'] 

## create a kmf object
kmf = KaplanMeierFitter() 

## Fit the data into the model
kmf.fit(durations, event_observed,label='Kaplan Meier Estimate')

## Create an estimate
kmf.plot(ci_show=False) ## ci_show is meant for Confidence interval, since our data set is too tiny, thus i am not showing it.

