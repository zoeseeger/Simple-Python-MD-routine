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
    d_no = 3                            # FOR X, Y AND Z
    p_no = len(coords[0])               # TO ITERATE THROUGH ATOMS
    
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
            mu, sigma = 0, 0.01          # CONSIDER SIGMA MORE CAREFULLY
            vel = np.random.normal(mu, sigma, size=(d_no, p_no))
            for row in vel:
                sum_velocity = 0
                for val in row:
                    sum_velocity = sum_velocity + val
                if sum_velocity > 0.2:  # CONSIDER CUTOFF MORE CAREFULLY, WHAT IS TOO LARGE?
                    print("SUM OF VELOCITIES ",sum_velocity, " TOO LARGE")

            # TEMPERATURE FROM VELOCITY MEAN SQUARED
            # NOT SURE ABOUT MATHS TO FIND active_temp
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
            # VELOCITIES USED TO FIND PREVIOUS POSITIONS
            pos_prev = np.zeros([d_no, p_no)]
            for i in d_no:
                for j in p_no:
                    pos_prev[i,j] = coords[i][j] - vel[i,j] * dt

            # EMPTY ACCELERATIONS
            acc = np.zeros([d_no, p_no])

        else:
            # AFTER FIRST STEP UPDATE POSITION, VELOCITY, ACCELERATION USING update()
            coords, vel, acc = update (p_no, d_no, coords, vel, force, acc, dt, step)
        
        # CALCULATE FORCE, POTENTIAL AND KINETIC ENERGIES USING compute()
        force, potential, kinetic = compute (p_no, d_no, coords, vel, step)

        # INITIAL ENERGY e0 CALCULATED
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


# --------------------------------------
# UPDATES POSITIONS, VELOCITIES AND ACCELERATIONS USING VELOCITY 
# VERLET ALGORITHM
# -------------------------------------
def update (p_no, d_no, coords, vel, f, acc, dt, step):
    if step < 2:
        print("BEGINNING UPDATE")    

    # UPDATE POSITIONS
    for i in range(0, d_no):
        for j in range(0, p_no):
            coords[i][j] = float(coords[i][j]) + vel[i,j] * dt + 0.5 * acc[i,j] * dt * dt
    
    # UPDATE VELOCITIES
    for i in range(0, d_no):
        for j in range(0, p_no):
            rmass = 1/coords[5][j]
            vel[i,j] = vel[i,j] + 0.5 * dt * (f[i,j] * rmass + acc[i,j])
    
    # UPDATE ACCELERATIONS
    for i in range(0, d_no):
        for j in range(0, p_no):
            rmass = 1/coords[5][j]
            acc[i,j] = f[i,j] * rmass

    return coords, vel, acc

# -------------------------------------
# FUNCTION COMPUTES FORCES AND ENERGIES
# -------------------------------------
def compute (p_no, d_no, coords, vel, step):
    import numpy as np
    from math import sin, sqrt
    
    if step < 2:
        print('BEGINNING COMPUTE')
    pi_d15 = 3.141592653589793
    
    # ZERO MATRICES
    force  = np.zeros([d_no, p_no])
    rij    = np.zeros(d_no) # SHAPE (3,) 

    potential = 0.0

    # COMPUTE ENERGY AND FORCES
    for i in range (0, p_no):
        # ALL PARTICLES j WITH PARTICLE i 
        for j in range(0, p_no):
            if (i != j): # IF i DOES NOT EQUAL j
                # COMPUTE Rij = DISPLACEMENT VECTOR - DIFFERENT IN X,Y,Z
                for k in range (0, d_no):
                    # ADDING TO rij [3,1] ONLY USED TO FIND d
                    rij[k] = float(coords[k][i]) - float(coords[k][j])
                
                # COMPUTE DISTANCE (d) AND TRUNCATED DISTANCE (d2)
                d = 0.0
                for k in range(0, d_no):
                    d = d + rij[k] ** 2
                # sqrt(d) = sqrt(X**2 + Y**2 + Z**2) = DISTANCE BETWEEN
                d = sqrt(d)
                d2 = min(d, pi_d15/2) # WILL RETURN WHICH EVER IS SMALLER
                # GIVE HALF OF THE POTENTIAL ENERGY TO PARTICLE j
                potential = potential + 0.5 * sin(d) * sin(d2)
                # ADD PARTICLE j'S CONTRIBUTION TO THE FORCE ON i
                for k in range(0, d_no):
                    # ADDING TO MATRIX force 
                    force[k,i] = force[k,i] - rij[k] * sin(2.0 * d2) / d

    # COMPUTE KINETIC ENERGY
    kinetic = 0.0
    for k in range(0, d_no):
        for j in range(0, p_no):
            # COORDS[5] = ATOMIC MASSES
            kinetic = kinetic + 0.5 * coords[5][j] * vel[k,j] ** 2
    
    return force, potential, kinetic

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





                



    

