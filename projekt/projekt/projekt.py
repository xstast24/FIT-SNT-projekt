from pylab import *
from scipy.integrate import solve_ivp


class Simulation:
    parameters = {
        'A_aorta': 4.15,
        'A_lv': 12,
        'b_index': [0, 435.2574, 4194.299, 3205.077, -1345.34],
        'BIS_0': 100,
        'BIS_MAX': 0,
        'c': 0.06263,
        'EC_50': 6.15e-5,
        'Eff_max_DP_E': 1.3,
        'Eff_max_DP_R': 0.5,
        'Eff_max_SNP_R': 0.635,
        'E_max0': 2.12,
        'f_R': 14.5,
        'g_index': [0, 24.456, 8.412, 4.667, 1.247],
        'k_20': 0.0093,
        'k_e0': 0.948,
        'K': 4.316,
        'MAP0': 90,
        'R_index': [1.59, 1.4, 2.92, 44.9, 44.9],
        'R_sys0': 0.0258,
        'V': 5000,
        'V_index': [2310, 7100, 11300, 3000, 5100],
        'V_lv': 85,
        'V_T': 500,  # dechovy objem [ml]
        'delta': 150,  # dead space [ml]
        'deltaQ': 300,  # ztraty [ml/min]
        'gama': 1.6,
        'ro': 1.05,
        'tau_DP': 2,
        'tau_SNP': 0.25,
        # dohledat
        # 'C_in': 1,#.01, stala davka 1 g/ml isoflurane
         # 'C_out': 0,
        # Q_in=p['f_R']*p['V_T']
        'Q_in': 7250,  # flow rate - vstupni objem (ml/min)
        # 5000->3500
        'Q_index': [3500, 2700, 510, 1100, 220]  # blood flow in compartment (ml/min) - hodnoty z http://anesthesiology.pubs.asahq.org/article.aspx?articleid=2035555
    }

    def run(self, t_span):
        """
        :param t_span: 2-tuple, interval of integration (t0, tf). The solver starts with t=t0 and integrates until it reaches t=tf.
        :return: Bunch object with the following fields defined: t, y, etc. -- see scipy/integrate/_ivp/ivp.py fo details
        """
        state0 = [0] * 17
        result = solve_ivp(self.do_step, t_span, state0, method='LSODA')  # , t_eval=np.linspace(0, 1.5, 15)) #method='LSODA',
        return result

    def do_step(self, t, state):
        PKM_state = state[0:16]
        PDM_state = state[16:]  # TODO
        PKM_results = self.pharmacokinetic_model(t, PKM_state, self.parameters)
        PDM_results = self.pharmacodynamic_model(t, PDM_state, self.parameters)
        return PKM_results + PDM_results

    def dosage(self, t):
        """nastaveni davkovani pro isofluran (C_in), dopamin (C_inf_DP) a nitroprusid sodny (C_inf_SNP) v danem case"""
        if (t > 5) and (t < 20):
            C_inf_DP = 1
        else:
            C_inf_DP = 0
        C_in = 0.01
        C_inf_SNP = 0
        return C_in, C_inf_DP, C_inf_SNP

    # pharmaco-kinetic model
    def pharmacokinetic_model(self, t, state, p):
        C_insp, C1_I, C2_I, C3_I, C4_I, C5_I, C1_DP, C2_DP, C3_DP, C4_DP, C5_DP, C1_SNP, C2_SNP, C3_SNP, C4_SNP, C5_SNP = state
        C_insp_old = C_insp
        C_index_I = [C1_I, C2_I, C3_I, C4_I, C5_I]
        C_index_DP = [C1_DP, C2_DP, C3_DP, C4_DP, C5_DP]
        C_index_SNP = [C1_SNP, C2_SNP, C3_SNP, C4_SNP, C5_SNP]

        C_in, C_inf_DP, C_inf_SNP = self.dosage(t)
        C_out = C_in - C_insp  # TODO

        # ISOFLURANE
        # respiratory system
        # koncentrace v plicich = (prisun isoflurance - naredeni - vyfouknuti) / objem plic
        C_insp = (p['Q_in'] * C_in - (p['Q_in'] - p['deltaQ']) * C_insp_old - p['f_R'] * (p['V_T'] - p['delta']) * (C_insp_old - C_out)) / p['V']

        # lungs - central compartment
        pom = 0
        for i in range(1, 5):
            pom += p['Q_index'][i] * (C_index_I[i] / p['R_index'][i] - C_index_I[0])
        # C1_I = (pom + p['f_R']*(p['V_T']-p['delta'])*(C_insp_old-C_index_I[0])) / p['V_index'][0]
        C1_I = (pom + p['f_R'] * (p['V_T'] - p['delta']) * (C_insp - C_index_I[0])) / p['V_index'][0]

        # liver - 2nd compartment
        C2_I = (p['Q_index'][1] * (C_index_I[0] - C_index_I[1] / p['R_index'][1]) -
                p['k_20'] * C_index_I[1] * p['V_index'][1]) / p['V_index'][1]

        # muscles, other organs and tissues, fat tissues - 3rd-5th compartment
        pom = [0, 0, 0, 0, 0]
        for i in range(2, 5):
            pom[i] = (p['Q_index'][i] * (C_index_I[0] - C_index_I[i] / p['R_index'][i])) / p['V_index'][i]
        C3_I = pom[2]
        C4_I = pom[3]
        C5_I = pom[4]

        # DOPAMINE
        # heart - central compartment
        pom = 0
        for i in range(1, 5):
            pom += p['Q_index'][i] * (C_index_DP[i] / p['R_index'][i] - C_index_DP[0])
        C1_DP = (pom + C_inf_DP - (C_index_DP[0] * p['V_index'][0]) / p['tau_DP']) / p['V_index'][0]

        # liver, muscles, other organs and tissues, fat tissues - 2nd-5th compartment
        pom = [0, 0, 0, 0, 0]
        for i in range(1, 5):
            pom[i] = (p['Q_index'][i] * (C_index_DP[0] - C_index_DP[i] / p['R_index'][i]) -
                      (C_index_DP[i] * p['V_index'][i]) / p['tau_DP']) / p['V_index'][i]
        C2_DP = pom[2]
        C3_DP = pom[2]
        C4_DP = pom[3]
        C5_DP = pom[4]

        # SODIUM NITROPRUSSIDE
        # heart - central compartment
        pom = 0
        for i in range(1, 5):
            pom += p['Q_index'][i] * (C_index_SNP[i] / p['R_index'][i] - C_index_SNP[0])
        C1_SNP = (pom + C_inf_SNP - (C_index_SNP[0] * p['V_index'][0]) / p['tau_SNP']) / p['V_index'][0]

        # liver, muscles, other organs and tissues, fat tissues - 2nd-5th compartment
        pom = [0, 0, 0, 0, 0]
        for i in range(1, 5):
            pom[i] = (p['Q_index'][i] * (C_index_SNP[0] - C_index_SNP[i] / p['R_index'][i]) -
                      (C_index_SNP[i] * p['V_index'][i]) / p['tau_SNP']) / p['V_index'][i]
        C2_SNP = pom[2]
        C3_SNP = pom[2]
        C4_SNP = pom[3]
        C5_SNP = pom[4]

        result = [C_insp, C1_I, C2_I, C3_I, C4_I, C5_I, C1_DP, C2_DP, C3_DP, C4_DP, C5_DP, C1_SNP, C2_SNP, C3_SNP, C4_SNP, C5_SNP]
        return result

    def pharmacodynamic_model(self, t, state, p):
        # todo
        return [0]


simulation = Simulation()
time_span = [0, 50]  # starting and final times
sol = simulation.run(time_span)

vykreslit_I = []
vykreslit_DP = []
vykreslit_SNP = []
for i in range(sol.t.size):
    prom = []
    # for j in range(int(sol.y.size/sol.t.size)):
    for j in range(6):
        prom.append(sol.y.item(j, i))
    vykreslit_I.append(prom)

    prom = []
    for j in range(6, 11):
        prom.append(sol.y.item(j, i))
    vykreslit_DP.append(prom)

    prom = []
    for j in range(11, 16):
        prom.append(sol.y.item(j, i))
    vykreslit_SNP.append(prom)

# efekt isofluranu na MAP
vykreslit_MAP = []
for i in range(sol.t.size):
    pom = 0
    for j in range(1, 5):
        pom += simulation.parameters['g_index'][j] * (1 + simulation.parameters['b_index'][j] * vykreslit_I[i][j + 1])
    vykreslit_MAP.append(simulation.parameters['Q_index'][0] / pom)

plt.figure(1)

plt.subplot(111)
plt.plot(sol.t, vykreslit_MAP)
plt.xlabel('Time (min)')
plt.ylabel('MAP (mmHg)')
plt.title('Mean arterial pressure (MAP)')

plt.figure(2)

plt.subplot(311)
plt.plot(sol.t, vykreslit_I)
plt.xlabel('Time (min)')
plt.ylabel('Concentration (g/mL)')
plt.title('Concentration of isoflurane')
plt.legend(('$C_{insp}$ (g/mL)', '$C_1$ (g/mL)', '$C_2$ (g/mL)', '$C_3$ (g/mL)', '$C_4$ (g/mL)', '$C_5$ (g/mL)'))

plt.subplot(312)
plt.plot(sol.t, vykreslit_DP)
plt.xlabel('Time (min)')
plt.ylabel('Concentration (g/mL)')
plt.title('Concentration of dopamine')
plt.legend(('$C_1$ (g/mL)', '$C_2$ (g/mL)', '$C_3$ (g/mL)', '$C_4$ (g/mL)', '$C_5$ (g/mL)'))

plt.subplot(313)
plt.plot(sol.t, vykreslit_SNP)
plt.xlabel('Time (min)')
plt.ylabel('Concentration (g/mL)')
plt.title('Concentration of sodium nitroprusside')
plt.legend(('$C_1$ (g/mL)', '$C_2$ (g/mL)', '$C_3$ (g/mL)', '$C_4$ (g/mL)', '$C_5$ (g/mL)'))

plt.show()
