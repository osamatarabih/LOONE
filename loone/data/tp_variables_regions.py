import os
import pandas as pd
from loone.utils import load_config


class TP_Variables:
    """Class representing TP variables."""
    def __init__(self, working_path: str):
        os.chdir(working_path)
        config = load_config(working_path)
        self.Z_sed = config['z_sed']
        self.Per_H2O_M = config['per_h2o_m']
        self.Per_H2O_S = config['per_h2o_s']
        self.Per_H2O_R = config['per_h2o_r']
        self.Per_H2O_P = config['per_h2o_p']
        self.N_Per = config['n_per']
        self.S_Per = config['s_per']
        ####
        self.Bulk_density_M = config['bulk_density_m']
        self.Bulk_density_S = config['bulk_density_s']
        self.Bulk_density_R = config['bulk_density_r']
        self.Bulk_density_P = config['bulk_density_p']
        ####
        self.Particle_density_M = config['particle_density_m']
        self.Particle_density_S = config['particle_density_s']
        self.Particle_density_R = config['particle_density_r']
        self.Particle_density_P = config['particle_density_p']
        ####
        self.A_Mud_N = config['a_mud_n']
        self.A_Mud_S = config['a_mud_s']
        self.A_Sand_N = config['a_sand_n']
        self.A_Sand_S = config['a_sand_s']
        self.A_Rock_N = config['a_rock_n']
        self.A_Rock_S = config['a_rock_s']
        self.A_Peat_N = config['a_peat_n']
        self.A_Peat_S = config['a_peat_s']
        self.A_N = self.A_Mud_N + self.A_Sand_N + self.A_Rock_N + self.A_Peat_N
        self.A_S = self.A_Mud_S + self.A_Sand_S + self.A_Rock_S + self.A_Peat_S
        A_tot = self.A_N + self.A_S
        self.Per_M_N = self.A_Mud_N / A_tot
        self.Per_M_S = self.A_Mud_S / A_tot
        self.Per_S_N = self.A_Sand_N / A_tot
        self.Per_S_S = self.A_Sand_S / A_tot
        self.Per_R_N = self.A_Rock_N / A_tot
        self.Per_R_S = self.A_Rock_S / A_tot
        self.Per_P_N = self.A_Peat_N / A_tot
        self.Per_P_S = self.A_Peat_S / A_tot

        self.Per_M_NN = self.A_Mud_N / self.A_N
        self.Per_M_SS = self.A_Mud_S / self.A_S
        self.Per_S_NN = self.A_Sand_N / self.A_N
        self.Per_S_SS = self.A_Sand_S / self.A_S
        self.Per_R_NN = self.A_Rock_N / self.A_N
        self.Per_R_SS = self.A_Rock_S / self.A_S
        self.Per_P_NN = self.A_Peat_N / self.A_N
        self.Per_P_SS = self.A_Peat_S / self.A_S
        self.Γ_inf = 91  # (mg/kg)
        #####Monthly
        self.v_burial_M = config["v_burial_m"]
        self.v_burial_S = config["v_burial_s"]
        self.v_burial_R = config["v_burial_r"]
        self.v_burial_P = config["v_burial_p"]

        # Read Calibration Outputs
        Cal_Res = pd.read_csv(config["nondominated_sol_var"])

        Par = Cal_Res["Par"]
        self.v_diff_M = Par[0]
        self.v_diff_S = Par[1]
        self.v_diff_R = Par[2]
        self.v_diff_P = Par[3]
        ####
        self.K_decomp_M = Par[4]
        self.K_decomp_S = Par[5]
        self.K_decomp_R = Par[6]
        self.K_decomp_P = Par[7]
        ###
        self.K_des_M = Par[8]
        self.K_des_S = Par[9]
        self.K_des_R = Par[10]
        self.K_des_P = Par[11]
        ####
        self.K_ads_M = Par[12]
        self.K_ads_S = Par[13]
        self.K_ads_R = Par[14]
        self.K_ads_P = Par[15]
        # v_settle = Par[16]
