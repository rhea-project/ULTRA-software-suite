## @ingroup Components-Costs
# Costs.py
#
# Created:
# Modified: Feb 2016, T. MacDonald
# Modified: Feb 2016, T. Orra

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------
from SUAVE.Components import Component
from SUAVE.Core import Data
from SUAVE.Methods.Costs.Correlations.Industrial_Costs import estimate_hourly_rates

# ----------------------------------------------------------------------
# Operating Costs class
# ----------------------------------------------------------------------
## @ingroup Components-Costs
class Operating_Costs(Data):
    """A class containing operating cost variables.
    
    Assumptions:
    None
    
    Source:
    N/A
    """    
    def __defaults__(self):
        """This sets the default values used in the operating cost methods.

        Assumptions:
        None

        Source:
        N/A

        Inputs:
        None

        Outputs:
        None

        Properties Used:
        N/A
        """          
        self.tag = 'operating_costs'
        self.depreciate_years = 0.0
        self.fuel_price       = 0.0
        self.oil_price        = 0.0
        self.insure_rate      = 0.0
        self.interest_rate    = 0.0
        self.maintenance_rate = 0.0
        self.pilot_rate       = 0.0
        self.crew_rate        = 0.0
        self.inflator         = 0.0
        self.reference_dollars= 0.0

        self.flight_cycles      = 0.0
        self.labor_costs        = 0.0
        self.repair_costs       = 0.0
        self.yearly_flight_time = 0.0
        self.navigation_fees    = 0.0

        self.electricity_price  = 0.0
        self.crew_complements   = 0.0

        # Capital costs
        self.capital= Data() 
        self.capital.spare_parts = Data()   
        self.capital.price       = Data()      
        self.capital.residual_value         = 0.0
        self.capital.residual_value_battery = 0.0

        self.capital.spare_parts.airframe       = 0.0
        self.capital.spare_parts.gas_turbine    = 0.0
        self.capital.spare_parts.electric_motor = 0.0
        self.capital.spare_parts.kS_PMAD        = 0.0

        self.capital.price.airframe         = 0.0
        self.capital.price.gas_turbine      = 0.0
        self.capital.price.electric_motor   = 0.0
        self.capital.price.PMAD             = 0.0
        self.capital.price.Battery          = 0.0
        self.capital.battery_cycles         = 0.0



# ----------------------------------------------------------------------
# Industrial Costs class
# ----------------------------------------------------------------------
## @ingroup Components-Costs
class Industrial_Costs(Data):
    """A class containing industrial cost variables.
    
    Assumptions:
    None
    
    Source:
    N/A
    """     
    def __defaults__(self):
        """This sets the default values used in the industrial cost methods.

        Assumptions:
        None

        Source:
        N/A

        Inputs:
        None

        Outputs:
        None

        Properties Used:
        N/A
        """  
        # inputs
        self.tag                            = 'industrial_costs'
        self.reference_year                 = 0.0
        self.production_total_units         = 0.0
        self.units_to_amortize              = None
        self.prototypes_units               = 0.0
        self.avionics_cost                  = 0.0
        self.test_facilities_cost           = 0.0
        self.manufacturing_facilities_cost  = 0.0
        self.escalation_factor              = 0.0
        self.development_total_years        = 0.0
        self.aircraft_type                  = None # ('military','general aviation','regional','commercial','business')
        self.difficulty_factor              = 1.0  # (1.0 for conventional, 1.5 for moderately advanc. tech., 2 for agressive use of adv. tech.)
        self.cad_factor                     = 1.0  # (1.2 for learning, 1.0 for manual, 0.8 for experienced)
        self.stealth                        = 0.0  # (0 for non-stealth, 1 for stealth)
        self.material_factor                = 1.0  # (1 for conventional Al, 1.5 for stainless steel, 2~2.5 for composites, 3 for carbon fiber)

        # hourly rates
        self.hourly_rates = Data()
        hourly_rates = self.hourly_rates
        hourly_rates.engineering     = 0.0
        hourly_rates.tooling         = 0.0
        hourly_rates.manufacturing   = 0.0
        hourly_rates.quality_control = 0.0

        return