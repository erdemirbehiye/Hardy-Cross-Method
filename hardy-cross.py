"""
#Pipe Diagram

   (0)           (2) 
A------------B--------C
|            |        | (5)
|        (3) |   (4)  |
|(1)         G--------H
|            |        | (9)
|            |(7)     |
D------------E--------F
   (6)           (8)


"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

def calculate_dQ(K,Q,indexes):
   
   total_dQ=[]
   
   for index in indexes:
      hf=np.empty(len(index))
      ccw=[]
      cw=[]
      for j in range(len(hf)):
         hf[j]=K[np.abs(index[j])]*pow(Q[np.abs(index[j])],2) #m
         if(index[j]<0):
            ccw.append([hf[j]])
         else:
            cw.append([hf[j]])
      
      np.seterr(divide='ignore', invalid='ignore')
      hf_Q=np.divide(hf,Q[np.abs(index)])#s/m^2
      hf_Q[np.isnan(hf_Q)] = 0.0
      hf_Q_ccw=[]
      hf_Q_cw=[]
      
      for k in range(len(hf_Q)):
         if index[k]<0:
            hf_Q_ccw.append([hf_Q[k]])
         else:
            hf_Q_cw.append([hf_Q[k]])
      
      numerator=np.sum(cw)-np.sum(ccw)
      denumerator=2*(np.sum(hf_Q_cw)+np.sum(hf_Q_ccw))
      dQ=numerator/denumerator

      total_dQ.append(dQ)
      prime_Q=[]
      
      if dQ<0:#ccw
         for l in range(len(index)):
            if index[l]<0:
               prime_Q.append(Q[abs(index[l])]-abs(dQ))
            else:
               prime_Q.append(Q[index[l]]+abs(dQ))
      else:#cw
         for l in range(len(index)):
            if index[l]<0:#-index
               prime_Q.append(Q[abs(index[l])]+abs(dQ))
            else:#+index
               prime_Q.append(Q[index[l]]-abs(dQ))

      for i,p in enumerate(abs(index)):
         Q[p] = prime_Q[i]

      for j in range(len(Q)):
         if Q[j]!=abs(Q[j]):
            for  k in indexes:
               if np.where(abs(k)==j):
                  a=np.where(abs(k)==j)
                  k[a]=-1*k[a]
                     
      Q=abs(Q)
      

   return total_dQ,Q,indexes

def check_criteria(dQ,criteria,K,Q,index):
   counter=0
   loop1=[]
   loop2=[]
   loop3=[]

   while abs(dQ)>=criteria:
      counter=counter+1
      dQs,Q,index=calculate_dQ(K,Q,index)
      loop1.append(abs(dQs[0]))
      loop2.append(abs(dQs[1]))
      loop3.append(abs(dQs[2]))

      check=[True for x in dQs if abs(x) < criteria]
      if len(check)==len(dQs):  
         print(dQs, ", this delta Q values are satisfied the criteria.")
         print("*"*50)
         break
         loop1=np.array(loop1)
         loop2=np.array(loop2)
         loop3=np.array(loop3)
   return Q,counter,loop1,loop2,loop3

def check_pressure(newQ,K,gama,Pa):
   p_criteria=185.0
   hf= K*pow(newQ,2)
   dp=Pa-hf[1]*gama-hf[6]*gama-hf[8]*gama
   if dp < p_criteria:
      print(dp," kPa is less than ",p_criteria,"  kPa, the pressure requirement at F won't be satisfied.")
   else:
      print(dp," kPa is greater than ",p_criteria,"  kPa, the pressure requirement at F is satisfied.")
   print("*"*60)
   return hf

f=np.array([0.019,0.020,0.021,0.021,0.021,0.021,0.021,0.022,0.021,0.022]) #friction factor
L=np.array([300,250,350,125,350,125,300,125,350,125]) #m
g=9.81
D=np.array([0.30,0.25,0.20,0.20,0.20,0.20,0.20,0.15,0.20,0.15])
e_D=np.array([0.00087,0.00104,0.00130,0.00130,0.00130,0.00130,0.00130,0.00173,0.00130,0.00173])
Q=np.array([0.2,0.1,0.08 ,0.12,0.02 ,0.03 ,0.1,0.0, 0.1,0.05])
K=[]

for i in range(len(f)):
   K.append(8*f[i]*L[i]/(pow(np.pi,2)*g*(pow(D[i],5))))#s^2/m^5

index_1=np.array([0,3,7,-1,-6])
index_2=np.array([2,5,-3,-4])
index_3=np.array([4,9,-7,-8])

index=[index_1,index_2,index_3]
f_index=[index_1,index_2,index_3] #indices for factorized parameters

dQ=1
criteria_1=0.005
criteria=pow(10,-6)
gama=9.790
Pa=gama*50

print("For Criteria is equal ",criteria_1)
_,counter,loop1,loop2,loop3=check_criteria(dQ,criteria_1,K,Q,index)
print(loop1)
print(loop2)
print(loop3)

counter=np.arange(counter)

"""#plot 
fig, ax = plt.subplots()
ax.yaxis.set_major_formatter(FormatStrFormatter('%.4f'))
ax.yaxis.set_ticks(np.arange(-2, 2, 0.001))

plt.axhline(y=0.005, color='r', linestyle='dashed')

plt.plot(counter,loop1)
plt.plot(counter,loop2)
plt.plot(counter,loop3)

plt.xlabel('Iteration') 
plt.ylabel('DQ Values Over Iteration') 
plt.legend(["Criteria ","DQ Values of Loop 1", "DQ Values of Loop 2","DQ Values of Loop 3"], loc ="upper right") 
#plt.savefig('hardy',dpi=150)

plt.show()
"""

print("For Criteria is equal ",criteria)
new_Q,counter,loop1,loop2,loop3=check_criteria(dQ,criteria,K,Q,index)
print(loop1)
print(loop2)
print(loop3)

print("Pressure Requirement at F")
hf=check_pressure(new_Q,K,gama,Pa)


"""counter=np.arange(counter)
#plot 
fig, ax = plt.subplots()
ax.yaxis.set_major_formatter(FormatStrFormatter('%.4f'))
ax.yaxis.set_ticks(np.arange(-2, 2, 0.001))

plt.axhline(y=1e-6, color='r', linestyle='dashed')

plt.plot(counter,loop1)
plt.plot(counter,loop2)
plt.plot(counter,loop3)

plt.xlabel('Iteration') 
plt.ylabel('DQ Values Over Iteration') 
plt.legend(["Criteria ","DQ Values of Loop 1", "DQ Values of Loop 2","DQ Values of Loop 3"], loc ="upper right") 
#plt.savefig('hardy',dpi=150)

plt.show()"""

#factorized parameters
Factor=20
f = np.array([x * Factor for x in f])
L = np.array([x * Factor for x in L])
D = np.array([x * Factor for x in D])
Q = np.array([x * Factor for x in Q])
p_criteria=185.0
K=[]
for i in range(len(f)):
   K.append(8*f[i]*L[i]/(pow(np.pi,2)*g*(pow(D[i],5))))

print("For Factorized parameters, the criteria is equal to ",criteria)
new_Q,c,loop1,loop2,loop3=check_criteria(dQ,criteria,K,Q,index)
print('DQs')
print(loop1)
print(loop2)
print(loop3)

c=np.arange(c)

#plot 
fig, ax = plt.subplots()
ax.yaxis.set_major_formatter(FormatStrFormatter('%.4f'))
ax.yaxis.set_ticks(np.arange(-2, 2, 0.01))

plt.axhline(y=1e-6, color='r', linestyle='dashed')

plt.plot(c,loop1)
plt.plot(c,loop2)
plt.plot(c,loop3)

plt.xlabel('Iteration') 
plt.ylabel('DQ Values Over Iteration') 
plt.legend(["Criteria ","DQ Values of Loop 1", "DQ Values of Loop 2","DQ Values of Loop 3"], loc ="upper right") 
#plt.savefig('hardy',dpi=150)

plt.show()

print("Pressure Requirement at F")
check_pressure(new_Q,K,gama,Pa)

head_Pump=10
delta_HG=(head_Pump-hf[4])*gama
p_F=Pa-hf[0]*gama-hf[3]*gama+delta_HG-hf[9]*gama


print("If pump added between the points G and H, the pressure requirement at F : ")
if p_F < p_criteria:
   print(p_F," kPa is less than ",p_criteria,"  kPa, the pressure requirement at F won't be satisfied.")
else:
   print(p_F," kPa is greater than ",p_criteria,"  kPa, the pressure requirement at F is satisfied.")

