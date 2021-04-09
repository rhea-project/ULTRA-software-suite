## @ingroup Input_Output-Results
# print_costs.py

# Modified: Aug 2020, S. Karpuk

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------
import SUAVE
import numpy as np

from SUAVE.Core import Units
import time                     # importing library
import datetime                 # importing library

# ----------------------------------------------------------------------
#  Methods
# ----------------------------------------------------------------------
## @ingroup Input_Output-Results
def print_costs(DOC,vehicle,filename='costs.dat'):
    """This creates a file showing aircraft costs

    Assumptions:
    None

    Source:
    N/A

    Inputs:

    Outputs:
    filename                  Saved file with name as above

    Properties Used:
    N/A
    """   

    # Unpack inputs
    flight_time = vehicle.mission_time
    range       = vehicle.range / Units.nautical_miles * 1.8        # in km
    Npax        = vehicle.passengers 

    DOC_Energy  = DOC.Energy 
    DOC_fuel    = DOC.Fuel_Energy 
    DOC_battery = DOC.Battery_Energy     
    DOC_Crew    = DOC.Crew        
    DOC_Ma      = DOC.Maintenance 
    DOC_Cap     = DOC.Capital     
    DOC_Fees    = DOC.Fees  
    DOC_tot     = DOC.Total      
    FC          = DOC.flight_cycles


    # write header of file
    fid = open(filename, 'w')  # Open output file
    fid.write('Output file with costs data \n\n')
    fid.write('  VEHICLE TAG          : ' + vehicle.tag + '\n\n') 
    fid.write('  DOC                  : ' + vehicle.tag + '\n') 
    fid.write('  Flight time          : ' + str(flight_time) + 'hr\n')
    fid.write('  Flight range         : ' + str(range) + 'km\n')
    fid.write('  Number of passengers : ' + str(Npax) + '\n\n')

    fid.write('  DOC break-down:\n')

    fid.write(   '  Parameter      | DOC/flight (EUR) | DOC/km (EUR)| DOC/hr (EUR) | DOC/hr/pax (EUR) | DOC/seat/mile (USD cent) \n')
    Parameter  = '  Fuel Energy    '                  + '|'   
    Energy     = str('%8.5f'   % DOC_fuel )           + '|'   
    DOC_ER     =  DOC_fuel /range
    Energy_km  = str('%8.5f'   % DOC_ER )             + '|' 
    DOC_hr     = DOC_fuel /range/flight_time
    Energy_hr  = str('%8.5f'   % DOC_hr)              + '|' 
    DOC_pax    = DOC_fuel /range/flight_time/Npax
    Energy_pax = str('%8.5f'   % DOC_pax)             + '|'
    DOC_seat_mile = DOC_ER / Npax * 1.6 / 0.72 * 100
    Energy_seat_mile = str('%8.5f'   % DOC_seat_mile) + '|'
    fid.write(Parameter + Energy + Energy_km + Energy_hr + Energy_pax + Energy_seat_mile + '\n')

    Parameter  = '  Battery Energy ' + '|'   
    Energy     = str('%8.5f'   % DOC_battery)     + '|'   
    DOC_ER     = DOC_battery/range
    Energy_km  = str('%8.5f'   % DOC_ER )     + '|' 
    DOC_hr     = DOC_battery/range/flight_time
    Energy_hr  = str('%8.5f'   % DOC_hr)     + '|' 
    DOC_pax    = DOC_battery/range/flight_time/Npax
    Energy_pax = str('%8.5f'   % DOC_pax)     + '|'
    DOC_seat_mile = DOC_ER / Npax * 1.6 / 0.72 * 100
    Energy_seat_mile = str('%8.5f'   % DOC_seat_mile)     + '|'
    fid.write(Parameter + Energy + Energy_km + Energy_hr + Energy_pax + Energy_seat_mile + '\n')

    Parameter  = '  Total Energy   '     + '|'   
    Energy     = str('%8.5f'   % DOC_Energy)     + '|'   
    DOC_ER     =  DOC_Energy/range
    Energy_km  = str('%8.5f'   % DOC_ER )     + '|' 
    DOC_hr     = DOC_Energy/range/flight_time
    Energy_hr  = str('%8.5f'   % DOC_hr)     + '|' 
    DOC_pax    = DOC_Energy/range/flight_time/Npax
    Energy_pax = str('%8.5f'   % DOC_pax)     + '|'
    DOC_seat_mile = DOC_ER / Npax * 1.6 / 0.72 * 100
    Energy_seat_mile = str('%8.5f'   % DOC_seat_mile)     + '|'
    fid.write(Parameter + Energy + Energy_km + Energy_hr + Energy_pax + Energy_seat_mile + '\n')

    Parameter  = '  Crew           '     + '|'   
    Crew     = str('%8.5f'   % DOC_Crew)    + '|'   
    Crew_ER  = DOC_Crew/range
    Crew_km  = str('%8.5f'   % Crew_ER)     + '|' 
    Crew_hr  = DOC_Crew/range/flight_time
    Crew_hr  = str('%8.5f'   % Crew_hr)     + '|' 
    Crew_pax = DOC_Crew/range/flight_time/Npax
    Crew_pax = str('%8.5f'   % Crew_pax)    + '|'
    DOC_seat_mile = Crew_ER / Npax * 1.6 / 0.72 * 100
    Crew_seat_mile = str('%8.5f'   % DOC_seat_mile)     + '|'
    fid.write(Parameter + Crew + Crew_km + Crew_hr + Crew_pax + Crew_seat_mile + '\n')

    Parameter  = '  Maintenance    '     + '|'   
    Ma      = str('%8.5f'   % DOC_Ma)     + '|' 
    Ma_ER   = DOC_Ma/range 
    Ma_km   = str('%8.5f'   % Ma_ER)     + '|' 
    Ma_hr   = DOC_Ma/range/flight_time
    Ma_hr   = str('%8.5f'   % Ma_hr)     + '|' 
    Ma_pax  = DOC_Ma/range/flight_time/Npax
    Ma_pax  = str('%8.5f'   % Ma_pax)     + '|'
    DOC_seat_mile = Ma_ER / Npax * 1.6 / 0.72 * 100
    Ma_seat_mile = str('%8.5f'   % DOC_seat_mile)     + '|'
    fid.write(Parameter + Ma + Ma_km + Ma_hr + Ma_pax + Ma_seat_mile + '\n')       

    Parameter  = '  Capital        '     + '|'   
    Cap     = str('%8.5f'   % DOC_Cap)     + '|'   
    Cap_ER  = DOC_Cap/range
    Cap_km  = str('%8.5f'   % Cap_ER)     + '|' 
    Cap_hr  = DOC_Cap/range/flight_time
    Cap_hr  = str('%8.5f'   % Cap_hr)     + '|' 
    Cap_pax = DOC_Cap/range/flight_time/Npax
    Cap_pax = str('%8.5f'   % Cap_pax)     + '|'
    DOC_seat_mile = Cap_ER / Npax * 1.6 / 0.72 * 100
    Cap_seat_mile = str('%8.5f'   % DOC_seat_mile)     + '|'
    fid.write(Parameter + Cap + Cap_km + Cap_hr + Cap_pax + Cap_seat_mile + '\n')   

    Parameter  = '  Fees           '     + '|'   
    Fees     = str('%8.5f'   % DOC_Fees)     + '|'   
    Fees_ER  = DOC_Fees/range
    Fees_km  = str('%8.5f'   % Fees_ER)     + '|' 
    Fees_hr  = DOC_Fees/range/flight_time
    Fees_hr  = str('%8.5f'   % Fees_hr)     + '|' 
    Fees_pax = DOC_Fees/range/flight_time/Npax
    Fees_pax = str('%8.5f'   % Fees_pax)     + '|'
    DOC_seat_mile = Fees_ER / Npax * 1.6 / 0.72 * 100
    Fees_seat_mile = str('%8.5f'   % DOC_seat_mile)     + '|'
    fid.write(Parameter + Fees + Fees_km + Fees_hr + Fees_pax + Fees_seat_mile + '\n')      

    Parameter  = '  Total          '     + '|'   
    tot     = str('%8.5f'   % DOC_tot)     + '|'   
    tot_ER  = DOC_tot/range
    tot_km  = str('%8.5f'   % tot_ER)     + '|' 
    tot_hr  = DOC_tot/range/flight_time
    tot_hr  = str('%8.5f'   % tot_hr)     + '|' 
    tot_pax = DOC_tot/range/flight_time/Npax
    tot_pax = str('%8.5f'   % tot_pax)     + '|'
    DOC_seat_mile = tot_ER / Npax * 1.6 / 0.72 * 100
    tot_seat_mile = str('%8.5f'   % DOC_seat_mile)     + '|'
    fid.write(Parameter + tot + tot_km + tot_hr + tot_pax + tot_seat_mile + '\n')  


    # close file
    fid.close  

    return

# ----------------------------------------------------------------------
#   Module Test
# ----------------------------------------------------------------------
if __name__ == '__main__':
    print(' Error: No test defined ! ')    
