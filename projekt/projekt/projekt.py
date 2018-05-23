from pylab import *
from scipy.integrate import solve_ivp
from enum import Enum
from math import sqrt


class Pos:
    EFF_DP_EMAX = 16
    EFF_DP_RSYS = 17
    EFF_SNP_RSYS = 18
    C_E2 = 19
    C_E3 = 20
    C_E4 = 21
    C_E5 = 22


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
        'k_1': 1,  # TODO vyzkouset, neni v clanku
        'k_2': 1,  # TODO vyzkouset, neni v clanku
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
        state0 = [0] * 20  # TODO
        return solve_ivp(self.do_step, t_span, state0, method='LSODA')  # , t_eval=np.linspace(0, 1.5, 15)) #method='LSODA',

    def do_step(self, t, state):
        PKM_results = self.pharmacokinetic_model(t, state, self.parameters)
        PDM_results = self.pharmacodynamic_model(t, state, self.parameters)
        # TODO jak toto funguje?
        #print(state)
        #print('\n')
        #print(PKM_results)
        #print(PDM_results)
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
        C_insp, C1_I, C2_I, C3_I, C4_I, C5_I, C1_DP, C2_DP, C3_DP, C4_DP, C5_DP, C1_SNP, C2_SNP, C3_SNP, C4_SNP, C5_SNP = state[0:16]
        C_insp_old = C_insp
        C_index_I = [C1_I, C2_I, C3_I, C4_I, C5_I]
        C_index_DP = [C1_DP, C2_DP, C3_DP, C4_DP, C5_DP]
        C_index_SNP = [C1_SNP, C2_SNP, C3_SNP, C4_SNP, C5_SNP]

        C_in, C_inf_DP, C_inf_SNP = self.dosage(t)
        C_out = C_in - C_insp  # TODO

        # ISOFLURANE
        # respiratory system
        # koncentrace v plicich = (prisun isoflurance - naredeni - vyfouknuti) / objem plic
        new_C_insp = (p['Q_in'] * C_in - (p['Q_in'] - p['deltaQ']) * C_insp_old - p['f_R'] * (p['V_T'] - p['delta']) * (C_insp_old - C_out)) / p['V']

        # lungs - central compartment
        pom = 0
        for i in range(1, 5):
            pom += p['Q_index'][i] * (C_index_I[i] / p['R_index'][i] - C_index_I[0])
        # C1_I = (pom + p['f_R']*(p['V_T']-p['delta'])*(C_insp_old-C_index_I[0])) / p['V_index'][0]
        new_C1_I = (pom + p['f_R'] * (p['V_T'] - p['delta']) * (new_C_insp - C_index_I[0])) / p['V_index'][0]

        # liver - 2nd compartment
        new_C2_I = (p['Q_index'][1] * (C_index_I[0] - C_index_I[1] / p['R_index'][1]) -
                p['k_20'] * C_index_I[1] * p['V_index'][1]) / p['V_index'][1]

        # muscles, other organs and tissues, fat tissues - 3rd-5th compartment
        pom = [0, 0, 0, 0, 0]
        for i in range(2, 5):
            pom[i] = (p['Q_index'][i] * (C_index_I[0] - C_index_I[i] / p['R_index'][i])) / p['V_index'][i]
        new_C3_I = pom[2]
        new_C4_I = pom[3]
        new_C5_I = pom[4]

        # DOPAMINE
        # heart - central compartment
        pom = 0
        for i in range(1, 5):
            pom += p['Q_index'][i] * (C_index_DP[i] / p['R_index'][i] - C_index_DP[0])
        new_C1_DP = (pom + C_inf_DP - (C_index_DP[0] * p['V_index'][0]) / p['tau_DP']) / p['V_index'][0]

        # liver, muscles, other organs and tissues, fat tissues - 2nd-5th compartment
        pom = [0, 0, 0, 0, 0]
        for i in range(1, 5):
            pom[i] = (p['Q_index'][i] * (C_index_DP[0] - C_index_DP[i] / p['R_index'][i]) -
                      (C_index_DP[i] * p['V_index'][i]) / p['tau_DP']) / p['V_index'][i]
        new_C2_DP = pom[2]
        new_C3_DP = pom[2]
        new_C4_DP = pom[3]
        new_C5_DP = pom[4]

        # SODIUM NITROPRUSSIDE
        # heart - central compartment
        pom = 0
        for i in range(1, 5):
            pom += p['Q_index'][i] * (C_index_SNP[i] / p['R_index'][i] - C_index_SNP[0])
        new_C1_SNP = (pom + C_inf_SNP - (C_index_SNP[0] * p['V_index'][0]) / p['tau_SNP']) / p['V_index'][0]

        # liver, muscles, other organs and tissues, fat tissues - 2nd-5th compartment
        pom = [0, 0, 0, 0, 0]
        for i in range(1, 5):
            pom[i] = (p['Q_index'][i] * (C_index_SNP[0] - C_index_SNP[i] / p['R_index'][i]) -
                      (C_index_SNP[i] * p['V_index'][i]) / p['tau_SNP']) / p['V_index'][i]
        new_C2_SNP = pom[2]
        new_C3_SNP = pom[2]
        new_C4_SNP = pom[3]
        new_C5_SNP = pom[4]

        result = [new_C_insp, new_C1_I, new_C2_I, new_C3_I, new_C4_I, new_C5_I, new_C1_DP, new_C2_DP, new_C3_DP, new_C4_DP, new_C5_DP, new_C1_SNP, new_C2_SNP, new_C3_SNP, new_C4_SNP, new_C5_SNP]
        return result

    def pharmacodynamic_model(self, t, state, p):
        C_insp, C1_I, C2_I, C3_I, C4_I, C5_I, C1_DP, C2_DP, C3_DP, C4_DP, C5_DP, C1_SNP, C2_SNP, C3_SNP, C4_SNP, C5_SNP = state[0:16]
        Eff_DP_Emax = state[Pos.EFF_DP_EMAX]
        Eff_DP_Rsys = state[Pos.EFF_DP_RSYS]
        Eff_SNP_Rsys = state[Pos.EFF_SNP_RSYS]
        # C_e2 = state[Pos.C_E2]
        # C_e3 = state[Pos.C_E3]
        # C_e4 = state[Pos.C_E4]
        # C_e5 = state[Pos.C_E5]

        # TODO je to C^N nebo C s hornim indexem N?
        new_Eff_DP_Emax = p['k_1'] * C1_DP * (p['Eff_max_DP_E'] - Eff_DP_Emax) - p['k_2'] * Eff_DP_Emax
        new_Eff_DP_Rsys = p['k_1'] * C1_DP * (p['Eff_max_DP_R'] - Eff_DP_Rsys) - p['k_2'] * Eff_DP_Rsys
        new_Eff_SNP_Rsys = p['k_1'] * C1_SNP * (p['Eff_max_SNP_R'] - Eff_SNP_Rsys) - p['k_2'] * Eff_SNP_Rsys
        E_max = p['E_max0'] * (1 + new_Eff_DP_Emax)
        R_sys = p['R_sys0'] * (1 - Eff_DP_Rsys - Eff_SNP_Rsys)

        d = (2*p['K']**2)**2 - 4 * (1 / R_sys**2) * (-2 * p['K']**2 * p['V_lv'] * E_max)  # diskriminant = b^2 - 4ac pro MAP rovnici
        if d < 0:
            MAP = 0
            print("Chyba - rovnice nema reseni")
        elif d == 0:
            MAP = (-(2*p['K']**2) + math.sqrt(d)) / (2 * (1 / R_sys**2))
        else:
            MAP = (-(2*p['K']**2) + math.sqrt(d)) / (2 * (1 / R_sys**2))  # koren = -b +- odmocnina z D to cele deleno 2a
            # druhy koren neni potreba: x2 = (-b - math.sqrt(d)) / (2 * a)

        # TODO BIS calculating - nejak mi neni jasne co tam presne znamenaji horni indexy gamma (degree of non-linearity). Rovnice nad nadpisem "2.3 baroreflex"
        # new_C_e2 = p['k_e0'] * (C1_I - C2_I)
        # new_C_e3 = p['k_e0'] * (C1_I - C3_I)
        # new_C_e4 = p['k_e0'] * (C1_I - C4_I)
        # new_C_e5 = p['k_e0'] * (C1_I - C5_I)

        return [new_Eff_DP_Emax, new_Eff_DP_Rsys, new_Eff_SNP_Rsys, MAP]


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

# MAP z pharmacokinetic modelu
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

# TODO jak je mozne, ze vykresluje linearne, kdyz MAP se meni?
plt.figure(3)
# MAP z pharmacodynamic modelu
plt.subplot(111)
plt.plot(sol.t, sol.y[19])
plt.xlabel('Time (min)')
plt.ylabel('MAP (mmHg)')
plt.title('TODO MAP - PDM')
plt.show()
