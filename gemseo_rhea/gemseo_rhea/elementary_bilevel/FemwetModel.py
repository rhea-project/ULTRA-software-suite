

class FemwetModel:

    def __init__(self):
        pass

    def execute(self, Z_planform, Z_stiffness, X_struct, X_mass_lump, X_aero, Y_dyn_load):

        stress_rf = 2*X_struct
        aileron_eff = 3*X_aero
        fuel_burnt = X_struct + X_aero + X_mass_lump + Z_planform + Z_stiffness + Y_dyn_load
        y_beam_model = 4*X_struct + X_mass_lump
        const_stiff = 5*X_struct

        return stress_rf, aileron_eff, fuel_burnt, y_beam_model, const_stiff


