import os
import pandas as pd
from loone.utils import load_config


class TP_Variables:
    """Class representing TP variables."""
    def __init__(self, working_path: str):
        os.chdir(working_path)
        config = load_config(working_path)

        # Sediment properties
        self.sediment_depth = config['z_sed']
        self.water_content_mud = config['per_h2o_m']
        self.water_content_sand = config['per_h2o_s']
        self.water_content_rock = config['per_h2o_r']
        self.water_content_peat = config['per_h2o_p']
        self.northern_percentage = config['n_per']
        self.southern_percentage = config['s_per']

        # Bulk density
        self.bulk_density_mud = config['bulk_density_m']
        self.bulk_density_sand = config['bulk_density_s']
        self.bulk_density_rock = config['bulk_density_r']
        self.bulk_density_peat = config['bulk_density_p']

        # # Particle density
        self.particle_density_mud = config['particle_density_m']
        self.particle_density_sand = config['particle_density_s']
        self.particle_density_rock = config['particle_density_r']
        self.particle_density_peat = config['particle_density_p']

        # # Area of sediments
        self.area_mud_north = config['a_mud_n']
        self.area_mud_south = config['a_mud_s']
        self.area_sand_north = config['a_sand_n']
        self.area_sand_south = config['a_sand_s']
        self.area_rock_north = config['a_rock_n']
        self.area_rock_south = config['a_rock_s']
        self.area_peat_north = config['a_peat_n']
        self.area_peat_south = config['a_peat_s']

        self.total_area_north = self.area_mud_north + self.area_sand_north + self.area_rock_north + self.area_peat_north
        self.total_area_south = self.area_mud_south + self.area_sand_south + self.area_rock_south + self.area_peat_south
        total_area = self.total_area_north + self.total_area_south

        # Percentage calculations
        self.percentage_mud_north = self.area_mud_north / total_area
        self.percentage_mud_south = self.area_mud_south / total_area
        self.percentage_sand_north = self.area_sand_north / total_area
        self.percentage_sand_south = self.area_sand_south / total_area
        self.percentage_rock_north = self.area_rock_north / total_area
        self.percentage_rock_south = self.area_rock_south / total_area
        self.percentage_peat_north = self.area_peat_north / total_area
        self.percentage_peat_south = self.area_peat_south / total_area

        self.percentage_mud_north_n = self.area_mud_north / self.total_area_north
        self.percentage_mud_south_s = self.area_mud_south / self.total_area_south
        self.percentage_sand_north_n = self.area_sand_north / self.total_area_north
        self.percentage_sand_south_s = self.area_sand_south / self.total_area_south
        self.percentage_rock_north_n = self.area_rock_north / self.total_area_north
        self.percentage_rock_south_s = self.area_rock_south / self.total_area_south
        self.percentage_peat_north_n = self.area_peat_north / self.total_area_north
        self.percentage_peat_south_s = self.area_peat_south / self.total_area_south

        self.inorganic_fraction = 91  # (mg/kg)

        # Burial velocities
        self.burial_velocity_mud = config["v_burial_m"]
        self.burial_velocity_sand = config["v_burial_s"]
        self.burial_velocity_rock = config["v_burial_r"]
        self.burial_velocity_peat = config["v_burial_p"]

        # Read Calibration Outputs
        calibration_results = pd.read_csv(config["nondominated_sol_var"])
        parameters = calibration_results["Par"]

        self.diffusion_velocity_mud = parameters[0]
        self.diffusion_velocity_sand = parameters[1]
        self.diffusion_velocity_rock = parameters[2]
        self.diffusion_velocity_peat = parameters[3]

        self.decomposition_rate_mud = parameters[4]
        self.decomposition_rate_sand = parameters[5]
        self.decomposition_rate_rock = parameters[6]
        self.decomposition_rate_peat = parameters[7]

        self.desorption_rate_mud = parameters[8]
        self.desorption_rate_sand = parameters[9]
        self.desorption_rate_rock = parameters[10]
        self.desorption_rate_peat = parameters[11]

        self.adsorption_rate_mud = parameters[12]
        self.adsorption_rate_sand = parameters[13]
        self.adsorption_rate_rock = parameters[14]
        self.adsorption_rate_peat = parameters[15]

        # self.settling_velocity = parameters[16]  # Uncomment if needed
