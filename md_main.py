#!/short/k96/zls565/installations/bin/python3
# USING PYTHON/3.4.3
# SEE: https://people.sc.fsu.edu/~jburkardt/py_src/md/md.py

# --------------------------
# PRODUCES ENERGIES AND PRINTS
# --------------------------
def md (coords, s_no = 500, dt = 0.1, temp = 300): # DEFAULT VALUES
    import numpy as np
    import sys

    sys = sys.version.split()
    d_no = 3                            # FOR X, Y AND Z, ROWS IN MATRIX
    p_no = len(coords[0])               # TO ITERATE THROUGH ATOMS, COLUMNS IN MATRIX
    
    # STARTING PARAMS
    print('')
    print('     Molecular Dynamics Simulator'                           )
#    print('     Written by Zoe Seeger, 2016'                            )
    print('     This programme is using Python version ' + sys[0]        )
    print('     With: '                                                 )
    print('         Spatial dimension   : ' , d_no                      )
    print('         Number of particles : ' , p_no                      )
    print('         Number of steps     : ' , s_no                      )
    print('         Time step (seconds) : ' , dt                        )
    print('')
    print('     Step        Potential           Kinetic         Relative Error')
    print('     --------------------------------------------------------------')

    step_print_index = 1
    step_print_num   = 100
    step_print       = 0 

    for step in range(0, s_no + 1):
        if (step == 0):
            # --------------
            # INITIATION
            # --------------            
            # RANDOM DISTR. VELOCITES ABOUT 0
            mu, sigma = 0, 0.01          # CONSIDER SIGMA MORE CAREFULLY ***
            vel = np.random.normal(mu, sigma, size=(d_no, p_no))
            for row in vel:
                sum_velocity = 0
                for val in row:
                    sum_velocity = sum_velocity + val
                if sum_velocity > 0.2:  # CONSIDER CUTOFF MORE CAREFULLY, WHAT IS TOO LARGE? ***
                    print("SUM OF VELOCITIES ",sum_velocity, " TOO LARGE")

            # TEMPERATURE FROM VELOCITY MEAN SQUARED
            # NOT SURE ABOUT MATHS TO FIND active_temp ***
            for i in d_no:
                for j in p_no:
                    active_temp = vel[i][j] ** 2 * coords[5][j] # T = m*v^2 / Nf * k_B
            # TIMESAVING?
            active_temp = active_temp / (3 * p_no * 1.38064852e-23) # Nf ~ 3N
            # SCALING FACTOR OF OLD TEMP -> NEW TEMP FOR VELOCITIES 
            scale_temp  = sqrt(temp/active_temp) 
            for i in d_no:
                for j in p_no:
                    vel[i,j] = vel[i,j] * scale_temp

            # EMPTY ACCELERATIONS
            acc = np.zeros([d_no, p_no])

        # UPDATE POSITION USING update_pos()
        coords = update_pos (p_no, d_no, coords, vel, force, acc, dt)
        
        # CALCULATE FORCE, POTENTIAL AND KINETIC ENERGIES USING compute()
        force = compute (p_no, d_no, coords, vel, dt)
        
        # UPDATE POSITION USING update_vel()
        vel, acc, kinetic = update_vel (p_no, d_no, coords, vel, force, acc, dt)
        
        # INITIAL ENERGY e0 CALCULATED FOR ERROR CALC
        if (step == 0):
            e0 = potential + kinetic

        # IF TRUE WILL PRINT STEP
        if (step == step_print):
            # CALCULATE RELATIVE ERROR WRT E0 AND PRINT STEP
            rel = (potential + kinetic - e0)/e0
            print('     %8d  %14f  %14f  %14g' % ( step, potential, kinetic, rel ))
            # INCREASE step_print_index
            step_print_index += step_print_index
            #print(step_print_num)
            # // IS INTEGER DIVISION 
            step_print = (step_print_index * s_no ) // step_print_num
        if (step == s_no):
            print('     %8d  %14f  %14f  %14g' % ( step, potential, kinetic, rel ))
    return

# -------------------------------------
# UPDATE POSITIONS
# -------------------------------------
def update_pos (p_no, d_no, coords, vel, acc, dt):
    if step < 2:
        print("BEGINNING UPDATE_POSITIONS")    

    # UPDATE POSITIONS
    for i in range(0, d_no):
        for j in range(0, p_no):
            # VELOCITY VERLET; TAYLOR SERIES EXPANSION FOR COORDS
            # r(t+dt) = r(t) + v(t)*dt + 1/2 * force(t)/mass * dt**2 
            # f(t) = a(t)/mass & f(t+dt) = force 
            coords[i][j] = float(coords[i][j]) + vel[i,j] * dt + 0.5 * acc[i,j] * dt * dt
    return coords

# -------------------------------------
# FUNCTION COMPUTES FORCES AND ENERGIES
# -------------------------------------
def compute (p_no, d_no, coords, vel, dt):
    import numpy as np
    from math import sin, sqrt
    
    if step < 2:
        print('BEGINNING COMPUTE')
    
    # ZERO MATRICES
    force  = np.zeros([d_no, p_no])

    # NEED TO CALL FUNCTION TO SUBMIT SPEC ***
    
    # READ IN FORCE DATA TO USE AS POTENTIAL ***
    
    # PULL FORCE/POTENTIAL ENERGY FOR ENERGY OF STEP ***
    # fx(r) = -dU(r)/dx        (partial derivative)
    # fx(r) = -(x/r)(dU(r)/dr) (partial derivative)
    
    return force
    
# --------------------------------------
# UPDATES VELOCITIES/ACC USING VELOCITY VERLET
# -------------------------------------
def update_vel (p_no, d_no, coords, vel, force, acc, dt):
    if step < 2:
        print("BEGINNING UPDATE_VEL")    
    
    # UPDATE VELOCITIES - VELOCITY VERLET TAYLOR SERIES EXPANSION FOR VELOCITIES
    # v(t+dt) = v(t) + (force(t) + force(t+dt)) * 1/2 * 1/mass * dt 
    for i in range(0, d_no):
        for j in range(0, p_no):
            rmass = 1/coords[5][j]
            vel[i,j] = vel[i,j] + 0.5 * dt * (f[i,j] * rmass + acc[i,j])
    
    # UPDATE ACCELERATIONS
    for i in range(0, d_no):
        for j in range(0, p_no):
            rmass = 1/coords[5][j]
            # f = ma
            acc[i,j] = f[i,j] * rmass
            
    
            
    # COMPUTE KINETIC ENERGY
    # USED TO FIND TOTAL ENERGY (+ POTENTIAL) AND FIND % ERROR IN CALC
    kinetic = 0.0
    for k in range(0, d_no):
        for j in range(0, p_no):
            # COORDS[5] = ATOMIC MASSES
            # 1/2 * m * v^2
            kinetic = kinetic + 0.5 * coords[5][j] * vel[k,j] ** 2

    return vel, acc, kinetic

# ------------------------
# TIME STAMP - NO PARAMS
# ------------------------
def t_stamp ():
    import time 
    t = time.time()
    print(time.ctime(t))
    return None

# -----------------------
# PRACTICE RUN OF MD SIMULATION
# -----------------------
def md_test ():
    from read_gms import xyz
    from time import clock 
    t_stamp()
    # COORDS: LIST OF LISTS; SYMBOL, CHARGE, X, Y, Z, MASS 
    coords = xyz()
    s_no = 100
    dt   = 0.1
    temp = 300
    wtime1 = clock()
    # THE MD SIMULATION 
    md(coords, s_no, dt, temp)
    wtime2 = clock()
    print('Simulation human time = ', wtime2 - wtime1)
    # TERMINATION
    print('')
    print('Normal Execution.')
    t_stamp()

md_test()





                



    

