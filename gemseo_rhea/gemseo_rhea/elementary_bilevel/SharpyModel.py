

class SharpyModel:

    def __init__(self):
        pass

    def execute(self, Z_planform, Y_beam_model):

        y_dyn_load = 2*Z_planform
        flutter = 3*Y_beam_model
        sum_mass_lump = Z_planform + Y_beam_model

        return y_dyn_load, flutter, sum_mass_lump