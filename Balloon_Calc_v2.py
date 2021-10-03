import math
import numpy as np
import pandas as pd
import random

e = 2.71828182846

def Height():
    H   =   list(range(0, 84000))
    return  H
#-----------------------------------------------
def Temperature(H):

    if H <= 11000:
        Tb      =   288.15
        Beta    =   -0.0065
        Hb      =   0
    elif H > 11000 and H <= 20000:
        Tb      =   216.65
        Beta    =   0
        Hb      =   11000
    elif H > 20000 and H <= 32000:
        Tb      =   216.65
        Beta    =   0.001
        Hb      =   20000
    elif H > 32000 and H <= 47000:
        Tb      =   228.65
        Beta    =   0.0028
        Hb      =   32000
    elif H > 47000 and H <= 51000:
        Tb      =   270.65
        Beta    =   0
        Hb      =   47000
    elif H > 51000 and H <= 71000:
        Tb      =   270.65
        Beta    =   -0.0028
        Hb      =   51000
    elif H > 71000 and H <= 84000:
        Tb      =   214.65
        Beta    =   -0.002
        Hb      =   71000
    Th   =   Tb + Beta * (H-Hb)

    return  Th
#-----------------------------------------------
def Pressure(H, Th):

    R = 287.052
    g = 9.80665

    if H <= 11000:
        Tb      =   288.15
        Beta    =   -0.0065
        Hb      =   0
        pb      =   101325
    elif H > 11000 and H <= 20000:
        Tb      =   216.65
        Beta    =   0
        Hb      =   11000
        pb      =   22632.1
    elif H > 20000 and H <= 32000:
        Tb      =   216.65
        Beta    =   0.001
        Hb      =   20000
        pb      =   5474.9
    elif H > 32000 and H <= 47000:
        Tb      =   228.65
        Beta    =   0.0028
        Hb      =   32000
        pb      =   868.02
    elif H > 47000 and H <= 51000:
        Tb      =   270.65
        Beta    =   0
        Hb      =   47000
        pb      =   110.91
    elif H > 51000 and H <= 71000:
        Tb      =   270.65
        Beta    =   -0.0028
        Hb      =   51000
        pb      =   66.939
    elif H > 71000 and H <= 84000:
        Tb      =   214.65
        Beta    =   -0.002
        Hb      =   71000
        pb      =   3.9564

    if Beta != 0:
        ph = pb * (1 + Beta/Tb * (H - Hb))**(-g/(Beta * R))

    elif Beta == 0:
        ph = pb * e**(-g/(R*Th) * (H - Hb))

    return ph
#-----------------------------------------------
def Density(Th, ph):
    R   =   287.052

    rho_h = ph / (R * Th)

    return rho_h
#-----------------------------------------------
def SpeedOfSound(Th):

    R   =   287.052
    k   =   1.41

    a_h = (k * R * Th)**(1/2)

    return a_h
#-----------------------------------------------
def MainLoop(wind_altitude_array, wind_direction_array, wind_speed_array):
    height  =   Height()
    GM      =   3.986004e+14
    ER      =   6.3781e+06
    SimBase =   0.1
    SimTime =   72000
    TimeVector = np.arange(0,SimTime, SimBase)
    T           =   [0] * len(TimeVector)
    p           =   [0] * len(TimeVector)
    rho         =   [0] * len(TimeVector)
    g           =   [0] * len(TimeVector)
    F_G         =   [0] * len(TimeVector)
    F_B         =   [0] * len(TimeVector)
    F_D         =   [0] * len(TimeVector)
    a           =   [0] * len(TimeVector)
    v           =   [0] * len(TimeVector)
    h           =   [0] * len(TimeVector)
    r           =   [0] * len(TimeVector)
    A_b         =   [0] * len(TimeVector)
    V_b         =   [0] * len(TimeVector)
    ascent_rate =   [0] * len(TimeVector)


    t_0 = 0
    rho_He = 0.1786
    C_d = 0.3
    #----------- INITIAL CONDITIONS ---------------
    V_b[0] = 1.1         # initial baloon volume
    BaloonMass = 0.45
    PayloadMass = 0.25
    r_pop = 2.36
    TotalMass = BaloonMass + PayloadMass
    initial_altitude = 100            # height above the mean sea level
    initial_latitude = 52.2
    initial_longitude = 51.1
    parachute_fall_rate_0 = 5
    #---------------------------------------------
    pop_time = 0
    pop_height = 0
    landing_time = 0
    flag = 0


    #pierwsza pętla
    h[0]                = 0 + initial_altitude
    T[0]                = Temperature(h[0])
    p[0]                = Pressure(h[0], T[0])
    rho[0]              = Density(T[0], p[0])
    k_constant          = (p[0]*V_b[0])/T[0]
    g[0]                = GM / (ER + h[0])**2
    F_G[0]              = TotalMass * g[0]
    F_D[0]              = 0
    V_b[0]              = (k_constant * T[0])/p[0]
    A_b[0]              = math.pi * ((3*V_b[0])/(4*math.pi))**(2/3)
    r[0]                = math.sqrt(A_b[0]/math.pi)
    F_B[0]              = -1 * V_b[0] * (rho_He - rho[0]) * g[0]
    a[0]                = 0
    v[0]                = 0
    calka_v             = 0
    calka_h             = 0 + h[0]
    latitude            = [0] * len(TimeVector)
    longitude           = [0] * len(TimeVector)
    north_velocity      = [0] * len(TimeVector)
    east_velocity       = [0] * len(TimeVector)
    wind_speed_test     = [0] * len(TimeVector)
    wind_direction_test = [0] * len(TimeVector)
    parachute_fall_rate = [0] * len(TimeVector)


    north_velocity[0]   = 0
    east_velocity[0]    = 0
    latitude[0]         = initial_latitude
    longitude[0]        = initial_longitude
    wind_speed          = 0  # m/s
    wind_direction      = 0  # deg

    for i in range(1,len(TimeVector)):
        T[i]            = Temperature(h[i-1])
        p[i]            = Pressure(h[i-1], T[i])
        rho[i]          = Density(T[i], p[i])
        V_b[i]          = (k_constant * T[i])/p[i]
        g[i]            = GM / (ER + h[i])**2

        if flag == 0:
            F_B[i]          = V_b[i] * rho[i] * g[i]
            F_G[i]          = TotalMass * g[i]
            A_b[i]          = math.pi * ((3*V_b[i])/(4*math.pi))**(2/3)
            r[i]            = math.sqrt(A_b[i]/math.pi)
            F_D[i]          = 0.5*C_d*A_b[i]*((v[i-1])**2)
            a[i]            = (F_B[i] - F_G[i] - F_D[i])/TotalMass
            calka_v         = calka_v + (a[i] + a[i-1]) * 0.5 * SimBase
            v[i]            = calka_v
            calka_h         = calka_h + (v[i] + v[i-1]) * 0.5 * SimBase
            h[i]            = calka_h


        if flag == 1:
            #v[i] = 8
            v[i] = parachute_fall_rate_0 * 0.4/ (math.sqrt(rho[i]/rho[0]))
            calka_h         = h[i - 1] - (v[i] + v[i-1]) * 0.5 * SimBase
            h[i]            = calka_h



        if wind_altitude_array[0] <= p[i] <= wind_altitude_array[1]:
            wind_speed     = wind_speed_array[0]
            wind_direction = wind_direction_array[0]

        elif wind_altitude_array[1] < p[i] <= wind_altitude_array[2]:
            wind_speed     = wind_speed_array[1]
            wind_direction = wind_direction_array[1]

        elif wind_altitude_array[2] < p[i] <= wind_altitude_array[3]:
            wind_speed     = wind_speed_array[2]
            wind_direction = wind_direction_array[2]

        elif wind_altitude_array[3] < p[i] <= wind_altitude_array[4]:
            wind_speed     = wind_speed_array[3]
            wind_direction = wind_direction_array[3]

        elif wind_altitude_array[4] < p[i] <= wind_altitude_array[5]:
            wind_speed     = wind_speed_array[4]
            wind_direction = wind_direction_array[4]

        elif wind_altitude_array[5] < p[i] <= wind_altitude_array[6]:
            wind_speed     = wind_speed_array[5]
            wind_direction = wind_direction_array[5]

        elif wind_altitude_array[6] < p[i] <= wind_altitude_array[7]:
            wind_speed     = wind_speed_array[6]
            wind_direction = wind_direction_array[6]

        elif wind_altitude_array[7] < p[i] <= wind_altitude_array[8]:
            wind_speed     = wind_speed_array[7]
            wind_direction = wind_direction_array[7]

        elif wind_altitude_array[8] < p[i] <= wind_altitude_array[9]:
            wind_speed     = wind_speed_array[8]
            wind_direction = wind_direction_array[8]

        elif wind_altitude_array[9] < p[i] <= wind_altitude_array[10]:
            wind_speed     = wind_speed_array[9]
            wind_direction = wind_direction_array[9]

        elif wind_altitude_array[10] < p[i] <= wind_altitude_array[11]:
            wind_speed     = wind_speed_array[10]
            wind_direction = wind_direction_array[10]

        elif wind_altitude_array[11] < p[i] <= wind_altitude_array[12]:
            wind_speed     = wind_speed_array[11]
            wind_direction = wind_direction_array[11]

        elif wind_altitude_array[12] < p[i] <= wind_altitude_array[13]:
            wind_speed     = wind_speed_array[12]
            wind_direction = wind_direction_array[12]

        elif wind_altitude_array[13] < p[i] <= wind_altitude_array[14]:
            wind_speed     = wind_speed_array[13]
            wind_direction = wind_direction_array[13]

        elif wind_altitude_array[14] < p[i] <= wind_altitude_array[15]:
            wind_speed     = wind_speed_array[14]
            wind_direction = wind_direction_array[14]

        elif wind_altitude_array[15] < p[i] <= wind_altitude_array[16]:
            wind_speed     = wind_speed_array[15]
            wind_direction = wind_direction_array[15]

        elif wind_altitude_array[16] < p[i] <= wind_altitude_array[17]:
            wind_speed     = wind_speed_array[16]
            wind_direction = wind_direction_array[16]

        elif wind_altitude_array[17] < p[i] <= wind_altitude_array[18]:
            wind_speed     = wind_speed_array[17]
            wind_direction = wind_direction_array[17]

        elif wind_altitude_array[18] < p[i] <= wind_altitude_array[19]:
            wind_speed     = wind_speed_array[18]
            wind_direction = wind_direction_array[18]

        elif wind_altitude_array[19] < p[i] <= wind_altitude_array[20]:
            wind_speed     = wind_speed_array[19]
            wind_direction = wind_direction_array[19]

        elif wind_altitude_array[20] < p[i] <= wind_altitude_array[21]:
            wind_speed     = wind_speed_array[20]
            wind_direction = wind_direction_array[20]

        elif wind_altitude_array[21] < p[i] < wind_altitude_array[22]:
            wind_speed     = wind_speed_array[21]
            wind_direction = wind_direction_array[21]

        elif p[i] >= wind_altitude_array[22]:
            wind_speed     = wind_speed_array[22]
            wind_direction = wind_direction_array[22]


        wind_speed_test[i] = wind_speed
        wind_direction_test[i] = wind_direction

        north_velocity[i]   = wind_speed * 0.514444444 * np.sin(3 * math.pi / 2 - math.radians(wind_direction))
        east_velocity[i]    = wind_speed * 0.514444444 * np.cos(3 * math.pi / 2 - math.radians(wind_direction))
        latitude[i]         = latitude[i - 1] + north_velocity[i] * SimBase * (360 / (2 * math.pi * (ER + h[i])))
        longitude[i]        = longitude[i - 1] + east_velocity[i] * SimBase * (360 / (2 * math.pi * (ER + h[i])))


        if r[i] >= r_pop and flag == 0:
            flag = 1
            pop_height  = h[i]
            calka_h = 0
            pop_time    = i * SimBase


        if h[i] <= initial_altitude and flag == 1:
            landing_time = i * SimBase
            break

    gravity_noise = [0] * len(TimeVector)
    buoyancy_noise = [0] * len(TimeVector)
    drag_noise = [0] * len(TimeVector)
    for i in range(0, i):
        gravity_noise[i] = F_G[i] + (random.random() - 0.5)
        buoyancy_noise[i] = F_B[i] + (random.random() - 0.5)
        drag_noise[i] = F_D[i] + (random.random() - 0.5)

    array1 = np.transpose(np.array([TimeVector, h, v, a, T, p, rho, g, F_G, F_B, F_D, V_b, A_b, r, latitude, longitude, gravity_noise, buoyancy_noise, drag_noise]))
    array1 = array1[0:i]
    array2 = array1[::200][:]

    output_data = pd.DataFrame(array2,columns=['Time [s]', 'Height [m]', 'Vertical speed [m/s]',
                                               'Acceleration [m/s^2]', 'Temperature [K]', 'Pressure [Pa]',
                                               'Density [kg/m^3]', 'Gravitational acceleration [m/s^2]',
                                               'Gravity [N]', 'Buoyancy [N]', 'Drag [N]',
                                               'Balloon volume [m^3]', 'Balloon cross-section area [m^2]', 'Baloon radius [m]',
                                               'Latitude [deg]', 'Longitude [deg]', 'Gravity Measured [N]', 'Buoyancy Measured [N]', 'Drag Measured [N]'])

    return output_data, pop_time, pop_height, north_velocity, east_velocity
#-----------------------------------------------

#-----------------------------------------------

def run():
    wind_forecast           = np.loadtxt('Wind_forecast_Central_EU_2.txt')
    wind_altitude_array     = wind_forecast[:, 0]
    wind_direction_array    = wind_forecast[:, 1]
    wind_speed_array        = wind_forecast[:, 2]

    for j in range(0, len(wind_altitude_array) - 1):
        wind_altitude_array[j]  = wind_altitude_array[len(wind_altitude_array) - 1 - j]
        wind_speed_array[j]     = wind_speed_array[len(wind_altitude_array) - 1 - j]
        wind_direction_array[j] = wind_direction_array[len(wind_altitude_array) - 1 - j]

    output_data, pop_time, pop_height, north_velocity, east_velocity = MainLoop(wind_altitude_array, wind_direction_array, wind_speed_array)
    return output_data, pop_time, pop_height, north_velocity, east_velocity


