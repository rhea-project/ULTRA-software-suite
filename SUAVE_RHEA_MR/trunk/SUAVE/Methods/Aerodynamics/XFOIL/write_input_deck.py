## @ingroup Methods-Aerodynamics-Xfoil
# write_input_deck.py
# 
# Created:  Aug 2019, S.Karpuk
# Modified: 


# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

from .purge_files import purge_files

## @ingroup Methods-Aerodynamics-AVL
def write_input_deck(xfoil_input,xfoil_object,M_eff,Cl_eff,Re_eff,flag,N_panels):
    """ This fucntions writes the execution steps used in the Xfoil call

    Assumptions:
        None
        
    Source:

    Inputs:
        xfoil_input
        xfoil_object
        M_eff                       [Unitless]
        Cl_eff                      [Unitless]
        Re_eff                      [Unitless]
        flag                        [String]
        N_panels                    [Unitless]

    Outputs:
        xfoil_input
 
    Properties Used:
        N/A
    """    

    # Set the filename up and Purge the file
    filename      = 'commands_' + xfoil_object.tag + '.dat'
    save_filename = 'polar_' + xfoil_object.tag + '.dat'
    purge_files([filename,save_filename])

    with open(filename, 'w') as deck_file:
        deck_file.write('PLOP\n')
        deck_file.write('G\n\n')
        deck_file.write('LOAD ' + xfoil_object.tag + '.dat\n')
        deck_file.write('PPAR\n')
        deck_file.write('N ' + str(N_panels) + '\n\n')
        deck_file.write('T 1\n\n\n')
        deck_file.write('OPER\n')
        deck_file.write('iter 200\n')
        deck_file.write('M ' + str(M_eff) + '\n')
        deck_file.write('VISC ' + str(round(Re_eff/1E5)*1E5) + '\n')
        deck_file.write('PACC\n')
        deck_file.write(save_filename + '\n\n')
        if flag == True:
            deck_file.write('CSEQ 0 ' + str(round(Cl_eff+0.1,1)) + ' 0.1\n\n\n')
        else:
            deck_file.write('Cl ' + str(Cl_eff) + '\n\n\n')
        deck_file.write('QUIT\n')

    xfoil_object.result_file  = save_filename
    xfoil_object.command_file = filename

    return xfoil_input
        
        
        
        
        
        
    
    
