import os
import subprocess
import itertools
import multiprocessing
import random

'''This script lunches a sweep of simulations exploring the behavior
of a disk packing.'''


def write_input_file(ac, seed, dirname):
    f = open(dirname + "/input_file", "w")
    f.write("""#nParticles           60
#seed                 """+seed+"""   (seed for random generator)
#box_w                10        (in disk units)
#box_h                10        (in disk units)
#freq                 80       (in Hz, for bottom movement)
#dimensionlessAc      """+ac+"""      (dimensionless acceleration of bottom)
#gravity              9.81     (m/s)
#gravityAngle         0.0      (in fractions of PI)
#bGamma               0        (bulk dissipation)
#timestep             1e-6     (in s, timestep for integrator)
#relaxTime            10       (in s, time for relaxation)
#thermalTime          100
#runTime              900      (in s, time for simulation)
#timeForGraph         1e0      (in s, time between graphics)
#timeForWriteRun      1e0      (in s, time between writes)
#timeForWriteThermal  1e9      (in s, time between writes)
#meanR                0.02     (in m, mean disk radius)
#density              3.57     (density of the material)
#kn                   4.5e6    (normal elastic constant)
#pr                   0.37     (poisson ratio of the material)
#mu                   0.1      (friction coefficient)
#vGamma               0.035    (viscoelastic dissipation)
#bCondType            1        (1: sinusoidal 2:random vib)
#relInitDisp          1        (between 1 and -1)
#wedge                0
""")


def run_simulation((ac, seed)):
    reRun = False
    execPath = "/home/moukarzel/ed-disks/disks-event"
    dirname = "seed."+seed+".ac."+ac
    if os.path.isdir(dirname):
        if not reRun:
            print("Skipping ac = ", ac, "with seed = ", seed,
                  ". Already been simulated.")
            return
    else:
        os.mkdir(dirname)
        # os.mkdir(dirname+"/Collisions")

    write_input_file(ac, seed, dirname)
    log = open(dirname + "/log", "w")
    errorLog = open(dirname + "/error_log", "w")
    # print("running ac = ", ac, "seed = ", seed)
    error = False
    if subprocess.call([execPath], cwd=dirname, stdout=log, stderr=errorLog):
        error = True
        print("Error at seed = ", seed, "ac = ", ac)
    log.close()
    errorLog.close()
    if not error:
        print("done ac = ", ac, "seed = ", seed)


def main():
    sweepdir = "ED-Sweep"
    if not os.path.isdir(sweepdir):
        os.mkdir(sweepdir)
    os.chdir(sweepdir)

    # List of acceleration amplitudes
    acList = ["00.5", "01.0", "02.0", "03.0", "04.0",
              "05.0", "06.0", "07.0", "08.0", "09.0", "10.0"]

    # List of seeds.
    random.seed(123456)
    seedList = [str(random.randint(100000, 900000)) for x in range(10)]

    pool_size = 16
    pool = multiprocessing.Pool(processes=pool_size)
    pool.map(run_simulation, itertools.product(acList, seedList))
    pool.close()  # no more tasks
    pool.join()  # wrap up current tasks
    os.chdir("..")

if __name__ == "__main__":
    main()
