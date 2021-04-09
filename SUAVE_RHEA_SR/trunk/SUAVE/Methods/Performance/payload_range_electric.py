## @ingroup Methods-Performance
# payload_range.py
#
# Created:  Apr 2014, T. Orra
# Modified: Jan 2016, E. Botero

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

from SUAVE.Core import Units, Data
import time
import numpy as np

# ----------------------------------------------------------------------
#  Calculate vehicle Payload Range Diagram
# ----------------------------------------------------------------------

## @ingroup Methods-Performance
def payload_range_electric(vehicle,mission,contingency,cruise_segment_tag,final_segment_tag,final_mission_tag,battery_TO,reserves=0.):
    """Calculates a vehicle's payload range diagram. Includes plotting.

    Assumptions:
    Constant altitude cruise

    Source:
    N/A

    Inputs:
    vehicle.mass_properties.
      operating_empty                     [kg]
      max_zero_fuel                       [kg]
      max_takeoff                         [kg]
      max_payload                         [kg]
      max_fuel                            [kg]
      takeoff                             [kg]
    mission.segments[0].analyses.weights.
      vehicle.mass_properties.takeoff     [kg]
    cruise_segment_tag                    <string>

    Outputs:
    payload_range.
      range                             [m]
      payload                           [kg]
      fuel                              [kg]
      takeoff_weight                    [kg]
    PayloadRangeDiagram.dat (text file)

    Properties Used:
    N/A
    """        
    # elapsed time start
    start_time = time.time()

    # Flags for printing results in command line, write output file, and plot
    iprint = 0      # Flag for print output data in the prompt line
    iwrite = 1      # Flag for write an output file
    iplot  = 1      # Flag for plot payload range diagram
    ### could be an user input.
    ##      output_type: 1: Print only              (light)
    ##      output_type: 2: Print + Write           (medium)
    ##      output_type: 3: Print + Write + Plot    (complete)

    #unpack
    masses     = vehicle.mass_properties
    if not masses.operating_empty:
        print("Error calculating Payload Range Diagram: Vehicle Operating Empty not defined")
        return True
    else:
        OEW = masses.operating_empty

    if not masses.max_zero_fuel:
        print("Error calculating Payload Range Diagram: Vehicle MZFW not defined")
        return True
    else:
        MZFW = vehicle.mass_properties.max_zero_fuel

    if not masses.max_takeoff:
        print("Error calculating Payload Range Diagram: Vehicle MTOW not defined")
        return True
    else:
        MTOW = vehicle.mass_properties.max_takeoff

    if not masses.max_payload:
        MaxPLD = MZFW - OEW  # If payload max not defined, calculate based in design weights
    else:
        MaxPLD = vehicle.mass_properties.max_payload
        MaxPLD = min(MaxPLD , MZFW - OEW) #limit in structural capability

    if not vehicle.propulsors.network.battery.mass_properties.mass: 
        print("Error calculating Payload Range Diagram: Vehicle Operating Empty not defined")
        return True
    else:
        Battery        = vehicle.propulsors.network.battery.mass_properties.mass  
        Battery_energy = vehicle.propulsors.network.battery.max_energy

    # Define payload range points
    #Point  = [ RANGE WITH MAX. PLD   , RANGE WITH MAX. FUEL , FERRY RANGE   ]
    TOW     = [ MTOW             , OEW              ]
    BATTERY = [ Battery_energy   , Battery_energy   ]
    PLD     = [ MaxPLD           , 0.               ]

    # allocating Range array
    R = [0,0]

    # evaluate the mission
    if iprint:
        print('\n\n\n .......... PAYLOAD RANGE DIAGRAM CALCULATION ..........\n')

     # Exclude reserves (works only for sizing scripts)
    vehicle.propulsors.network.battery.mass_properties.mass +=  reserves
    vehicle.store_diff()                                   # store the differnece if any

    # loop for each point of Payload Range Diagram
    for i in range(len(TOW)):
##    for i in [2]:
        if iprint:
            print(('   EVALUATING POINT : ' + str(i+1)))

        # Define takeoff weight
        mission.segments[0].analyses.weights.vehicle.mass_properties.takeoff = TOW[i]

        # Evaluate mission with current TOW
        results     = mission.evaluate()
        segment     = results.segments[cruise_segment_tag]

        # Distance convergency in order to have total fuel equal to target fuel
        #
        # User don't have the option of run a mission for a given fuel. So, we
        # have to iterate distance in order to have total fuel equal to target fuel
        #

        smgm_keys = len(results.segments.keys())

        maxIter = 10   # maximum iteration limit
        tol = 1     # fuel convsergency tolerance
        err = 9999.    # error to be minimized
        iter = 0       # iteration count

        while abs(err) > tol and iter < maxIter:
            iter = iter + 1

            # Difference between used energy and target energy
            missingCharge = results.segments[final_mission_tag].conditions.propulsion.battery_energy[-1,0] -  \
                            reserves * vehicle.propulsors.network.battery.specific_energy #+ battery_TO 
                            
            # Current distance and fuel consuption in the cruise segment
            CruiseDist   = np.diff( segment.conditions.frames.inertial.position_vector[[0,-1],0] )[0]        # Distance [m]
            CruiseEnergy = segment.conditions.propulsion.battery_energy[0,0] - segment.conditions.propulsion.battery_energy[-1,0] # [Ws]

            # Current specific range (m/Ws)
            CruiseSR    = CruiseDist / CruiseEnergy        # [m/Ws]

            # Estimated distance that will result in total energy used = target energy
            DeltaDist  =  CruiseSR *  missingCharge
            mission.segments[cruise_segment_tag].distance = (CruiseDist + DeltaDist)

            # running mission with new distance
            results = mission.evaluate()
            segment = results.segments[cruise_segment_tag]

            # Difference between burned fuel and target fuel
            err = results.segments[final_mission_tag].conditions.propulsion.battery_energy[-1,0] -  \
                  reserves * vehicle.propulsors.network.battery.specific_energy #+ battery_TO 

            if iprint:
                print(('     iter: ' +str('%2g' % iter) + ' (kg) | Residual : '+str('%8.0F' % err)))

        # Allocating resulting range in ouput array.
        R[i] = ( results.segments[final_segment_tag].conditions.frames.inertial.position_vector[-1,0] ) * Units.m / Units.km      #Distance [km]

    # Inserting point (0,0) in output arrays
    R.insert(0,0)
    PLD.insert(0,MaxPLD)
    TOW.insert(0,0)
    BATTERY.insert(-1,0)

    # packing results
    payload_range          = Data()
    payload_range.range     = R
    payload_range.payload        = PLD
    payload_range.takeoff_weight = TOW
    payload_range.reserves       = reserves

    # Write output file
    if iwrite:
        import datetime                 # importing library

        fid = open('PayloadRangeDiagram.dat','w')   # Open output file
        fid.write('Output file with Payload Range Diagram details\n\n') #Start output printing

        fid.write( ' Maximum Takeoff Weight ...........( MTOW ).....: ' + str( '%8.0F'   %   MTOW   ) + ' kg\n' )
        fid.write( ' Operational Empty Weight .........( OEW  ).....: ' + str( '%8.0F'   %   OEW    ) + ' kg\n' )
        fid.write( ' Maximum Zero Fuel Weight .........( MZFW ).....: ' + str( '%8.0F'   %   MZFW   ) + ' kg\n' )
        fid.write( ' Maximum Payload Weight ...........( PLDMX  )...: ' + str( '%8.0F'   %   MaxPLD ) + ' kg\n' )
        fid.write( ' Maximum Fuel Weight ..............( FUELMX )...: ' + str( '%8.0F'   %   Battery) + ' kg\n' )
        fid.write( ' Reserve Fuel  .................................: ' + str( '%8.0F'   %   reserves)+ ' kg\n\n' )

        fid.write( '    RANGE    |   PAYLOAD   |   BATTERY   |    TOW      |  \n')
        fid.write( '     km      |     kg      |    kg       |     kg      |  \n')

        for i in range(len(TOW)):
            fid.write( str('%10.0f' % R[i]) + '   |' + str('%10.0f' % PLD[i]) + '   |' + str('%10.0f' % BATTERY[i]) + '   |' + ('%10.0f' % TOW[i]) + '   |\n')

        # Print timestamp
        fid.write(2*'\n'+ 43*'-'+ '\n' + datetime.datetime.now().strftime(" %A, %d. %B %Y %I:%M:%S %p"))
        fid.close

    # Print data in command line
    if iprint:
        print( '\n\n                        RESULTS\n')
        print( '    RANGE    |   PAYLOAD   |   BATTERY   |    TOW      |')
        print( '     km      |     kg      |    kg       |     kg      |')
        for i in range(len(TOW)):
            print(( str('%10.0f' % R[i]) + '   |' + str('%10.0f' % PLD[i]) + '   |' + str('%10.0f' % BATTERY[i]) + '   |' + ('%10.0f' % TOW[i]) + '   |'))
        print(('\n\n   Elapsed time: ' + str('%6.2f' % (time.time() - start_time)) + 's'))

    #   Plot Payload Range
    if iplot:

        #import pylab
        import pylab as plt

        title = "Payload Range Diagram"
        plt.figure(0)
        plt.plot(R,PLD,'r')
        plt.xlabel('Range (km)'); plt.ylabel('Payload (kg)'); plt.title(title)
        plt.grid(True)
        plt.savefig("payload_range.png")
        #plt.show()

    return payload_range
