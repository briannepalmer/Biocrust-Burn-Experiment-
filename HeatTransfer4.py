# constants

bd = 1.4  # bulk density (g/cm3)
som = 0.07  # soil organic matter content (g/g)
vwct0 = 0.07  # initial volumetric water content (cm3/cm3)
gamma = 0.005  # thermal conductivity, W/cm*K
dt = 0.1  # time step, s
dz = 0.1  # soil layer thickness, cm
h = 2266.7  # latent heat of evaporation, J/g H2O (or cm3 H2O)
basaltemp = 20  # initial temperature (degrees C)
peaktemp = 500  # peak temperature
peakduration = 480  # time (s) at peak temperature
rampup = 240  # time to reach peak temperature
decay = 1800  # time to drop from peak back to basal
runtime = peakduration + rampup + decay
tfile = "temperature.txt"  # file names for model output
vwcfile = "vwc.txt"  # remember to change these for each run, or they will overwrite previous files


def surfacetemp(t):
    """returns surface temperature for a given time point"""
    from math import log, exp
    slope = (peaktemp - basaltemp) / rampup
    k = log(basaltemp / peaktemp) / decay
    if t <= rampup:
        return basaltemp + t * slope
    elif t < rampup + peakduration:
        return peaktemp
    else:
        return peaktemp * exp((t - rampup - peakduration) * k)


def spheat(vwc, om=som, p=bd):
    """specific heat as function of water, organic matter and bulk density"""
    minvol = p * (1 - om) / 2.65  # mineral soil by volume
    orgvol = p * om / 1.5  # om by volume
    return 4.184 * (0.45 * minvol + 0.6 * orgvol + vwc)  # J/cm3*K


# initialize matrices
ts0 = [surfacetemp(0)]  # soil temp matrix, 0,1 = previous and current time steps
ts1 = [surfacetemp(dt)]
vwc0 = [0]  # soil water content matrix, 0,1 = previous and current time steps
vwc1 = [0]
for i in range(200):
    ts0.append(basaltemp)
    ts1.append(basaltemp)
    vwc0.append(vwct0)
    vwc1.append(vwct0)

time = 0
counter = 0
tf = open(tfile, 'w')
vwcf = open(vwcfile, 'w')
tf.write('time,')
vwcf.write('time,')
print('time')
for col in range(201):  # print column headers for output files (depth in cm)
    tf.write(str(round(col * dz, 1)))
    tf.write(',')
    vwcf.write(str(round(col * dz, 1)))
    vwcf.write(',')
    print(str(round(col * dz, 1)))
tf.write('\n')
vwcf.write('\n')
print('\n')
while time < runtime:
    if counter % 100 == 0:  # every 100 time points, write to screen and output files
        print(round(time, 1))
        tf.write(str(round(time, 1)))
        tf.write(',')
        for elem in ts0:
            print('{elem:.1f}')  # I'm learning how to format - round looks way easier
            tf.write(str(round(elem, 1)))
            tf.write(',')
        print("\nVWC")
        tf.write('\n')
        vwcf.write(str(round(time, 1)))
        vwcf.write(',')
        for wc in vwc0:
            print(round(wc, 2))
            vwcf.write(str(round(wc, 2)))
            vwcf.write(',')
        print("\n")
        vwcf.write('\n')
        counter = 0
    for i in range(199):
        if ts0[i + 1] < 100 or vwc0[i + 1] < 0.000001:  # vwc started rising again due to discrete discontinuities
            hevap = 0  # heat of evaporation term
            vwc1[i + 1] = vwc0[i + 1]
            ch = spheat(vwc1[i + 1])
        else:
            hflux = gamma * ((ts0[i + 2] - ts0[i + 1]) - (ts0[i + 1] - ts0[i])) / dz  # calculate heat flux
            vwc1[i + 1] -= dt * hflux / h  # convert to water loss
            if vwc1[i + 1] < 0:  # (don't go negative)
                vwc1[i + 1] = 0
            ch = spheat(vwc1[i + 1])
            hevap = (h / ch) * (vwc1[i + 1] - vwc0[i + 1]) / dt  # calculate hevap from dvwc/dt
        ts1[i + 1] = ts0[i + 1] + dt * (gamma / ch) * ((ts0[i + 2] - ts0[i + 1]) - (ts0[i + 1] - ts0[i])) / dz + hevap
    time += dt
    counter += 1
    ts1[0] = surfacetemp(time)
    ts1[200] = basaltemp  # defined by boundary condition, increase model depth if necessary
    ts0 = ts1.copy()
    vwc0 = vwc1.copy()
tf.close()
vwcf.close()
