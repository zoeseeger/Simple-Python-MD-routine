#!/short/k96/zls565/installations/bin/python3
# USING PYTHON/3.4.3
# SOURCE: https://people.sc.fsu.edu/~jburkardt/py_src/md/md.py
# SOURCE: Frenkel, Smit, Understanding Molecular Simulation: From Algorithms to Applications

# --------------------------
# PRODUCES ENERGIES AND PRINTS
# --------------------------
def md (coords, s_no, dt, temp, original_f): # DEFAULT VALUES
    import numpy as np
    import math
    import os
    import sys

    sys  = sys.version.split()
    d_no = 3                            # FOR X, Y AND Z, ROWS IN MATRIX
    print(coords[0][0])
    p_no = len(coords[0])               # TO ITERATE THROUGH ATOMS, COLUMNS IN MATRIX

    # STARTING PARAMS
    print(''                                                                   )
    print('     Molecular Dynamics Simulator'                                  )
    print('     Written by Zoe Seeger, 2016'                                   )
    print('     This programme is using Python version ' + sys[0]              )
    print('     With: '                                                        )
    print('         Spatial dimension   : ' , d_no                             )
    print('         Number of particles : ' , p_no                             )
    print('         Number of steps     : ' , s_no                             )
    print('         Time step (seconds) : ' , dt                               )
    print(''                                                                   )
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
            active_temp = 0
            for i in range(0, d_no):
                for j in range(0, p_no):
                    # T = m*v^2 / Nf * k_B ~ IS THIS CORRECT?
                    active_temp = active_temp + vel[i,j] ** 2 * coords[5][j]
            # TIMESAVING?
            active_temp = active_temp / (3 * p_no * 1.38064852e-23) # Nf ~ 3N
            # SCALING FACTOR OF OLD TEMP -> NEW TEMP FOR VELOCITIES
            scale_temp  = math.sqrt(temp/active_temp)
            for i in range(0, d_no):
                for j in range(0, p_no):
                    vel[i,j] = vel[i,j] * scale_temp

            # EMPTY ACCELERATIONS
            acc = np.zeros([d_no, p_no])

            # TO STORE MD SIMS
            os.mkdir("gamess_sub")

        # UPDATE POSITION USING update_pos()
        coords = update_pos (p_no, d_no, coords, vel, acc, dt, step)

        # CALCULATE FORCE, POTENTIAL AND KINETIC ENERGIES USING compute()
        force, potential = compute (p_no, d_no, coords, vel, dt, step, original_f)

        # UPDATE POSITION USING update_vel()
        vel, acc, kinetic = update_vel (p_no, d_no, coords, vel, force, acc, dt, step)

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
def update_pos (p_no, d_no, coords, vel, acc, dt, step):
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
def compute (p_no, d_no, coords, vel, dt, step, original_f):
    import numpy as np
    import os
    import re

    if step < 2:
        print('BEGINNING COMPUTE')

    # ZERO MATRICES
    force  = np.zeros([d_no, p_no])

    # SUBMIT SPEC
    os.mkdir("gamess_sub/step_" + str(step))

    # SAVE ALL TEXT EXCEPT COORDS IN xyz
    save_coords = True
    inp         = []
    ut          = []
    for line in original_f:
        if save_coords:
            inp.append(line)
        if re.search('END', line):
            save_coords = True
        elif not save_coords:
            line_mod = line.split()
            ut.append([line_mod[0], line_mod[1]])
        if re.search('FMOXYZ', line):
            save_coords = False
    # ADD NEW COORDS TO xyz
    for i in range(0, d_no):
        for j in range(0, p_no):
            ut[j].append(float(coords[i][j]))
    print(ut)

    lines_job = [
                 "#!/bin/bash"                   ,
                 "#PBS -P k96 "                  ,
                 "#PBS -l walltime=4:00:00"      ,
                 "#PBS -l ncpus=16"              ,
                 "#PBS -l mem=40GB"              ,
                 "#PBS -l jobfs=60GB"            ,
                 "#PBS -l wd"                    ,
                 ""                              ,
                 "module unload openmpi/1.6.3"   ,
                 "module load openmpi/1.8.4"     ,
                 "/short/k96/apps/gamess/rungms.rika water.inp $PBS_NCPUS 01 > water.log"
                                                  ]

    with open("gamess_sub/step_" + str(step) + "/.inp", 'w+') as f:
        f.writelines(inp)
        for line in ut:
            for val in line:
                f.write(" " + str(val) + "\t")
            f.write("\n")
        f.write(" $END")
    with open("gamess_sub/step_" + str(step) + "/.job", 'w+') as f:
        f.writelines(lines_job)


    # READ IN FORCE DATA TO USE AS POTENTIAL ***
    potential = 0

    # PULL FORCE/POTENTIAL ENERGY FOR ENERGY OF STEP ***
    # fx(r) = -dU(r)/dx        (partial derivative)
    # fx(r) = -(x/r)(dU(r)/dr) (partial derivative)

    return force, potential

# --------------------------------------
# UPDATES VELOCITIES/ACC USING VELOCITY VERLET
# -------------------------------------
def update_vel (p_no, d_no, coords, vel, force, acc, dt, step):
    if step < 2:
        print("BEGINNING UPDATE_VELOCITIES")

    # UPDATE VELOCITIES - VELOCITY VERLET TAYLOR SERIES EXPANSION FOR VELOCITIES
    # v(t+dt) = v(t) + (force(t) + force(t+dt)) * 1/2 * 1/mass * dt
    for i in range(0, d_no):
        for j in range(0, p_no):
            rmass = 1/float(coords[5][j])
            vel[i,j] = vel[i,j] + 0.5 * dt * (force[i,j] * rmass + acc[i,j])

    # UPDATE ACCELERATIONS
    for i in range(0, d_no):
        for j in range(0, p_no):
            rmass = 1/float(coords[5][j])
            # f = ma
            acc[i,j] = force[i,j] * rmass


    # COMPUTE KINETIC ENERGY
    # USED TO FIND TOTAL ENERGY (+ POTENTIAL) AND FIND % ERROR IN CALC
    kinetic = 0.0
    for k in range(0, d_no):
        for j in range(0, p_no):
            # COORDS[4] = ATOMIC MASSES
            # 1/2 * m * v^2
            kinetic = kinetic + 0.5 * float(coords[5][j]) * vel[k,j] ** 2

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
    # COORDS: LIST OF LISTS; X, Y, Z, SYMBOL, CHARGE, MASS
    original_f, coords = xyz()
    s_no = 100
    dt   = 0.1
    temp = 300
    wtime1 = clock()
    # THE MD SIMULATION
    md(coords,  s_no, dt, temp, original_f)
    wtime2 = clock()
    print('Simulation human time = ', wtime2 - wtime1)
    # TERMINATION
    print('')
    print('Normal Execution.')
    t_stamp()

md_test()
