# READ IN GAMESS INPUT FILE

def xyz():
    import sys
    import re
    from periodictable import C,H,Cl,B,F,N 
    import numpy as np
    
    # USE "script.py filename.inp"
    #file = sys.argv[1]
    file = "c1mim-cl-p10.inp"
    
    # COORDS SAVED 
    save_coords = False
    coords      = []
    # COORDS IS LIST OF LISTS; [ATM SYMBOL, ATOMIC NUMBER, X, Y, Z, APPEND MASS]
    with open(file, 'r+') as f:
        for line in f:
            if re.search('END', line):
                save_coords = False
            elif save_coords:
                line = line.split()
                coords.append([line[2], line[3],line[4],line[0],line[1]])    
            elif re.search('FMOXYZ', line):
                save_coords = True
    # COORDS TRANSPOSED 
    coords = list(map(list, zip(*coords)))
    # CHARGE == ATOMIC NUMBER,
    # MAKE SECOND LETTER OR MORE LOWERCASE SO CAN BE READ BY periodictable
    coords.append([])
    for symbol in coords[3]:
        symbol2 = symbol
        if len(symbol2) > 1:
            for lett in symbol:
                if lett == symbol[0]:
                    symbol = lett 
                else:
                    symbol = symbol + lett.lower()
        # APPEND MOLAR MASS
        coords[5].append(eval(symbol).mass)

    return coords


