#constants
global G,m_sun_kg, c, sbc,pi, h,k 
G=6.67430e-11 #'''Nm2/kg2'''
m_sun_kg = 1.988e30 #'''kg'''#defined again in functioning
c=299792458 #'''m/s'''
sbc=5.67e-8 #watt/m2K4
pi=3.141592653589
h=6.62607015e-34 #J/Hz
k=1.380649e-23 #m2 kg /(Ks2)

FO1=3.2928e15
FO2=8.4922e15

# r_{s,ab} in units of R_S
R_AB_LOOKUP = {
    (1e8, 0.01): 50.54,
    (1e8, 0.05): 192.53,
    (5e8, 0.01): 59.92,
    (5e8, 0.05): 226.28,
    (1e9, 0.01): 64.51,
    (1e9, 0.05): 242.41,}
