from pylab import *
from scipy.integrate import solve_ivp

parameters={
'A_aorta': 4.15,
'A_lv': 12,
'b_index': [0,435.2574,4194.299,3205.077,-1345.34],
'BIS_0': 100,
'BIS_MAX': 0,
'c': 0.06263,
'EC_50': 6.15e-5,
'Eff_max_DP_E': 1.3,
'Eff_max_DP_R': 0.5,
'Eff_max_SNP_R': 0.635,
'E_max0': 2.12,
'f_R': 14.5,
'g_index': [0,24.456,8.412,4.667,1.247],
'k_20': 0.0093,
'k_e0': 0.948,
'K': 4.316,
'MAP0': 90,
'R_index':[1.59,1.4,2.92,44.9,44.9],
'R_sys0': 0.0258,
'V': 5000,
'V_index': [2310,7100,11300,3000,5100],
'V_lv': 85,
'V_T': 500, #dechovy objem [ml]
'delta': 150, #dead space [ml]
'deltaQ': 300, #ztraty [ml/min]
'gama': 1.6,
'ro': 1.05,
'tau_DP': 2,
'tau_SNP': 0.25,

#dohledat
'C_in': 1,#.01, stala davka 1 g/ml isoflurane
'C_out': 0,
# Q_in=p['f_R']*p['V_T']
'Q_in': 7250,  #flow rate - vstupni objem (ml/min)

#5000->3500
'Q_index': [5000,2700,510,1100,220]#blood flow in compartment (ml/min) - hodnoty z http://anesthesiology.pubs.asahq.org/article.aspx?articleid=2035555
}

#def Pharmacokinetic(t,state,p):
def Pharmacokinetic(t,state,p):
  C_insp,C1_I,C2_I,C3_I,C4_I,C5_I=state
  C_insp_old=C_insp
  C_index_I=list([C1_I,C2_I,C3_I,C4_I,C5_I])

  #ISOFLURANE
  #respiratory system
  p['C_out']=p['C_in']-C_insp #TODO
  #koncentrace v plicich = (prisun isoflurance - naredeni - vyfouknuti) / objem plic
  C_insp=(p['Q_in']*p['C_in'] - (p['Q_in'] - p['deltaQ'])*C_insp_old - p['f_R']*(p['V_T'] - p['delta'])*(C_insp_old - p['C_out']))/p['V']
  
  #lungs - central compartment
  pom=0
  for i in range(1,5):
    pom+=p['Q_index'][i]*(C_index_I[i]/p['R_index'][i]-C_index_I[0])
  #C1_I = (pom + p['f_R']*(p['V_T']-p['delta'])*(C_insp_old-C_index_I[0])) / p['V_index'][0]
  C1_I = (pom + p['f_R']*(p['V_T']-p['delta'])*(C_insp_old-C_index_I[0])) / p['V_index'][0]

  #liver - 2nd compartment
  C2_I = (p['Q_index'][1]*(C_index_I[0]-C_index_I[1]/p['R_index'][1]) -
          p['k_20']*C_index_I[1]*p['V_index'][1])/p['V_index'][1]
  
  #muscles, other organs and tissues, fat tissues - 3rd-5th compartment
  pom=[0,0,0,0,0]
  for i in range(2,5):
    pom[i]=(p['Q_index'][i]*(C_index_I[0]-C_index_I[i]/p['R_index'][i]))/p['V_index'][i]
  C3_I=pom[2]
  C4_I=pom[3]
  C5_I=pom[4]

  result=[C_insp,C1_I,C2_I,C3_I,C4_I,C5_I]
  return result


#t = np.linspace(0, 10, 101)

state0 = [0,0,0,0,0,0]

sol = solve_ivp(lambda t,y: Pharmacokinetic(t,y,parameters), [0,60], state0) #, t_eval=np.linspace(0, 1.5, 15)) #method='LSODA', 

vykreslit=[]
for i in range(sol.t.size):
  prom=[]
  for j in range(int(sol.y.size/sol.t.size)):
  #for j in [0,1]:
    prom.append(sol.y.item(j,i))
  vykreslit.append(prom)

plt.plot(sol.t, vykreslit)
xlabel('TIME (min)')
ylabel('Concentration (g/mL)')
title('Concentration of isoflurane inspired by the patient')
legend(('$C_{insp}$ (g/mL)','$C_1$ (g/mL)','$C_2$ (g/mL)','$C_3$ (g/mL)','$C_4$ (g/mL)','$C_5$ (g/mL)'))
plt.show()