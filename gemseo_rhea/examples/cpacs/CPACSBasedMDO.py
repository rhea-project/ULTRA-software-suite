
import numpy as np

from gemseo_rhea.cpacs import MDODisciplineCPACSBased


class TestDiscipline(MDODisciplineCPACSBased):

    def _run(self):
        x1 = self.get_inputs_by_name("x1")
        y1 = np.sin(x1) - np.exp(x1)
        self.local_data["y1"] = y1

    def _compute_jacobian(self, inputs=None, outputs=None):
        self._init_jacobian(with_zeros=True)
        x1 = self.get_inputs_by_name("x1")
        y = np.cos(x1) - np.exp(x1)
        self.jac["y1"] = {}
        self.jac["y1"]["x1"] = np.atleast_2d(y)


if __name__ == "__main__":

    from gemseo_rhea.cpacs import CPACSDesignSpace , CPACSStructureData

    # build CPACS data structure and select input variables
    input_cpacs_file = "input_wing.xml"
    cpacs_data = CPACSStructureData(input_cpacs_file)

    # select from generic XPath...
    cpacs_data.select_input_from_xpath(xpath="./vehicles/aircraft/model/wings/wing[@uID='wing1']/"
                                             +"positionings/positioning[@uID='pos3']/sweepAngle",
                                       name="x1")
    # ... or select from dedicated function
    # cpacs_data.select_sweep_angle_variable(name="x1",
    #                                        wing_ID="wing1",
    #                                        position_ID="pos3")

    # add a cpacs output
    cpacs_data.select_output_from_xpath("./vehicles/aircraft/model/analyses/aeroelastics/myResponse", "y1")

    # build design space from existing mapping
    design_space = CPACSDesignSpace(cpacs_data.input_mapping)

    # define variable data
    design_space.set_variable_data("x1", l_b=-2., u_b=2., value=-1.)

    # build discipline
    test_disc = TestDiscipline(cpacs_data.input_mapping, cpacs_data.output_mapping)

    from gemseo.api import create_scenario
    sc = create_scenario(test_disc, "DisciplinaryOpt", "y1", design_space, scenario_type="MDO")

    sc.execute({"algo":"SLSQP", "max_iter":100})
    res = sc.get_optimum()
    print(res)

    test_disc.set_cpacs_values()
    cpacs_data.write_xml("output_wing_discipline.xml")
