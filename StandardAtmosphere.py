import math
import numpy as np

e   =   2.71828182846
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
    elif i > 71000 and i <= 84000:
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
def AscentLoop(height):
    GM      =   3.986004e+14
    ER      =   6.3781e+06
    SimBase =   0.1
    SimTime =   84
    TimeVector = np.arange(0,SimTime, SimBase)
    g           =   [0] * len(TimeVector)
    F_G         =   [0] * len(TimeVector)
    F_B         =   [0] * len(TimeVector)
    F_D         =   [0] * len(TimeVector)
    A_b         =   [1] * len(TimeVector)  # change
    V_b         =   [0] * len(TimeVector)
    ascent_rate =   [0] * len(TimeVector)
    
    burst_altitude = 0
    time_to_burst = 0
    neck_lift = 0
    launch_radius = 0
    launch_volume = 0
    t_0 = 0
    mass_b = 0
    mass_p = 0
    rho_He = 0.1786
    C_d = 0.3
    BaloonMass = 1
    PayloadMass = 1
    TotalMass = BaloonMass + PayloadMass

    #burst_volume = (4.0/3.0) * math.pi * pow(balloon_burst_diameter / 2.0, 3)

    for i in range(len(TimeVector)):
        g[i]            = GM / (ER + height[i])**2
        F_B[i]          = V_b[i] * (rho_He - Pressure(height[i], Temperature(height[i]))) * g[i]
        F_G[i]          = TotalMass * g[i]
        #V_b[i]          = 
        #A_b[i]          = 
        #ascent_rate[i]   = math.sqrt(2*(F_B[i] - F_G[i])/(C_d * Pressure(height[i], Temperature(height[i])) * A_b[i]))
        #F_D[i]           = C_d * Pressure(height[i], Temperature(height[i])) * pow(ascent_rate,2) * A_b[i] / 2

    
                                
    return F_G, TimeVector
#-----------------------------------------------

#-----------------------------------------------
height          =   Height()

F_G, TimeVector = AscentLoop(height)
print(F_G)
print(TimeVector)
