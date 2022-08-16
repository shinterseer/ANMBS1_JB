import numpy as np

def p_Sat_calc(temp):
    '''
    Magnus formula to calculate the saturation vapourpressure based on the temperature.
    :param temp: Temperature in °C
    :return: p_Sat in Pa
    '''
    return 0.61094 * np.exp((17.625 * temp) / (temp + 243.04)) * 1000


def phi_kelvin_calc(p_suc):
    '''
    Kelvin formula to calculate the relative humidity based on the suction pressure in the pores.
    :param p_suc: (in Pa?????)
    :return: relative humidity (dimensionless)
    '''
    rho_w = 1000  # kg/m3
    Rv = 462  # J/ (kgK) Spezifische Gaskonstante für Wasserdampf
    T_ref = 293.15  # K
    rh = np.exp(-p_suc / (rho_w * Rv * T_ref))
    return rh


def p_suc_kelvin_calc(rh):
    '''
    Inverse of Kelvin formula for phi
    :param phi: 0.0 ... 1.0
    :return: p_suc (in Pa?????)
    '''

    rho_w = 1000  # kg/m3
    Rv = 462  # J/ (kgK) Spezifische Gaskonstante für Wasserdampf
    T_ref = 293.15  # K

    p_suc = - rho_w * Rv * T_ref * np.log(rh)
    return p_suc


def avg_conductivity_harmonic(c1, c2):
    return 2 * c1 * c2 / (c1 + c2)

def avg_conductivity_linear(c1, c2):
    return .5 * (c1 + c2)

def get_material(mat_string):
    materials = {
        'AAC_A4': {'description': 'Konstanten von Andi nach Bednar Habil A4',
                     'free_saturation': 350.,
                     'mu': 6.,
                     'pore_size': 1.0e-6,
                     'n': 5.,
                     'A': 6.2
                    },
        'AAC_A4_mod': {'description': 'Von Simon modifizierte Konstanten von Andi nach Bednar Habil A4',
                     'free_saturation': 350.,
                     'mu': 6.,
                     'pore_size': 1.0e-6,
                     'n': 7,
                     'A': 6.
                    },
        'AAC_A4_mod_dry': {'description': 'Von Simon modifizierte Konstanten von Andi nach Bednar Habil A4 - modifiziert für trocknung - siehe Caption zu Figure 9',
                     'free_saturation': 350.,
                     'mu': 4.,
                     'pore_size': 1.0e-6,
                     'n': 2.,
                     'A': 6.
                    },
        'Ziegel1': {'description': 'Ziegel Falkenlöwe lt. masea-ensan.com',
                     'free_saturation': 239.,
                     'mu': 10.7,
                     'pore_size': 3.0e-6, # hab ich optisch hingefittet aus sorptionsisotherme
                     'n': 13., # ich glaub, das sollen wir nehmen?
                     'A': 15.6
                    },
        'Porenbeton': {'description': 'Porenbeton YTONG Nord lt. masea-ensan.com',
                    'free_saturation': 320.,
                    'mu': 8.87,
                    'pore_size': 3.0e-6,  # übernommen von Ziegel1
                    'n': 13.,  # ich glaub, das sollen wir nehmen?
                    'A': 2.35
                    }
    }
    return materials[mat_string]


class one_d_problem:

    def __init__(self, res, sim_time, material, length, init_w=0, w_west=1, w_east=0):
        '''
        res ... resolution of spacial discretization (number of cells)
        sim_time ... simulated time in hours
        material ... string; must match an entry in get_material()
        length ... length of the specimen in meters
        init_w ... initial value: initial water content 0: dry, 1: freely saturated
        w_west ... boundary condition: water content of left border 0: dry, 1: freely saturated
        w_east ... boundary condition: water content of right border 0: dry, 1: freely saturated
        '''

        # Simulation Parameters
        # ---------------------
        self.resolution = res + 2
        self.sim_time = sim_time  # Simulation time in hours
        self.fluid_conduction = True  # should fluid water conduction be simulated
        self.vapour_diffusion = True # should vapour diffusion be simulated

        # Discretization Scheme
        # -----------------------
        self.length = length  # sample size in meters
        self.dx = self.length / (self.resolution - 2)
        self.x = np.linspace(self.dx / 2, self.length - self.dx / 2, num=self.resolution - 2)

        # Material
        # ----------------------------
        temp_mat = get_material(material)
        self.material = temp_mat['description']
        self.free_saturation = temp_mat['free_saturation']
        self.A = temp_mat['A']  # kg/m**2h**0.5  Moisture uptake coefficient
        self.mu = temp_mat['mu']  # Wasserdampf - Diffusionwiderstandszahl der Probe bei phi = 0
        self.n = temp_mat['n']
        self.pore_size = temp_mat['pore_size']

        self.delta_p_0 = 1.5 * 10 ** (-6)  # kg/mhPa
        self.delta_p = self.delta_p_0 / self.mu

        temperature = 23
        self.p_sat = np.full(self.resolution, p_Sat_calc(temperature), dtype=np.float64)

        # Initial Conditions
        # --------------------
        self.w_zero = 10 ** -9  # Numerical zero for water content
        self.w = np.full(self.resolution, max(init_w * self.free_saturation, self.w_zero), dtype=np.float64)

        #        BC
        # --------------
        # default problem: all flows enabled
        self.fluid_flow_west = True
        self.fluid_flow_east = True
        self.vapour_flow_west = True
        self.vapour_flow_east = True

        self.w[0] = max(w_west * self.free_saturation, self.w_zero)
        self.w[-1] = max(w_east * self.free_saturation, self.w_zero)


    def w_calc(self, P_suc):
        '''
        Function to calculate the water content based on the suction pressure in the cell.
        :param P_suc:
        :return: Water content of cell
        '''
        return self.free_saturation / (1. + self.pore_size * P_suc)


    def K_w_calc(self, p_suc):
        dwdp_suc = (-1) * self.free_saturation * self.pore_size / (1 + self.pore_size * p_suc)**2
        K_w = - dwdp_suc * (self.n + 1.) / (2. * self.n) * (self.A / self.free_saturation)**2 * (
              self.w_calc(p_suc) / self.free_saturation) ** self.n * (
              self.n + 1 - (self.w_calc(p_suc) / self.free_saturation)**self.n)
        return K_w


    def P_suc_calc(self, w):
        '''
        Inverse of w_calc().
        :param w: Water content of cell.
        :param free_saturation: Free Saturation of material
        :param pore_size: Free fitting parameter. Pore size of material.
        :return:
        '''
        P_suc = (self.free_saturation / w - 1.) / self.pore_size

        return P_suc


    # Function to solve - this is w_dot = f(t, w)
    # ----------------------------------------------
    def dwdt_calc(self, t, w):
        # Numerical approximation to avoid divsion by zero
        w[w == 0] = self.w_zero
        # Recalculate all dependant variables
        p_suc = self.P_suc_calc(w)
        K_w = self.K_w_calc(p_suc)
        rh = phi_kelvin_calc(p_suc)
        p_vap = rh * self.p_sat[0] # this means, temperature is assumed to be constant
        dwdt = np.zeros(w.shape)

        # define index-ranges for readability
        # -----------------------------------
        W = range(1, w.shape[0]-3) # west
        P = range(2, w.shape[0]-2) # point
        E = range(3, w.shape[0]-1) # east

        # middle Part
        # ------------
        # liquid water conduction
        K_W = avg_conductivity_linear(K_w[W], K_w[P])
        K_E = avg_conductivity_linear(K_w[P], K_w[E])
        dwdt[P] += 1/ self.dx**2 * ( (K_W + K_E) * p_suc[P] - K_W*p_suc[W] - K_E*p_suc[E]) * self.liquid_conduction

        # vapour diffusion
        dwdt[P] += (-1) * self.delta_p / self.dx**2 * (2*p_vap[P] - p_vap[W] - p_vap[E]) * self.vapour_diffusion

        # Western Boundary
        # ------------------
        # liquid water conduction
        K_W1 = avg_conductivity_linear(K_w[0], K_w[1])
        K_E1 = avg_conductivity_linear(K_w[1], K_w[2])
        dwdt[1] += (-1) * 1/ self.dx**2 * K_W1 * (p_suc[0] - p_suc[1]) * self.fluid_flow_west * self.liquid_conduction
        dwdt[1] += 1/ self.dx**2 * K_E1 * (p_suc[1] - p_suc[2]) * self.liquid_conduction

        # vapour diffusion
        dwdt[1] += (-1) * self.delta_p / self.dx**2 * (p_vap[1] - p_vap[0]) * self.vapour_flow_west * self.vapour_diffusion
        dwdt[1] += self.delta_p / self.dx**2 * (p_vap[2] - p_vap[1]) * self.vapour_diffusion

        # Eastern Boundary
        # ------------------
        # liquid water conduction
        K_Wend = avg_conductivity_linear(K_w[-3], K_w[-2])
        K_Eend = avg_conductivity_linear(K_w[-2], K_w[-1])
        dwdt[-2] += (-1) * 1/ self.dx**2 * K_Wend * (p_suc[-3] - p_suc[-2]) * self.liquid_conduction
        dwdt[-2] += 1/ self.dx**2 * K_Eend * (p_suc[-2] - p_suc[-1]) * self.fluid_flow_east * self.liquid_conduction

        # vapour diffusion
        dwdt[-2] += (-1) * self.delta_p / self.dx**2 * (p_vap[-2] - p_vap[-3]) * self.vapour_diffusion
        dwdt[-2] += self.delta_p / self.dx**2 * (p_vap[-1] - p_vap[-2]) * self.vapour_flow_east * self.vapour_diffusion

        return dwdt
