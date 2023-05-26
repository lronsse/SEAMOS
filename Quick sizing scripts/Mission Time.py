# Total Mission time for 10 units for the Puffin, Hybrid and Multi-system

# aerial distance to BWZ III
BWZ_dist = 30000  # m
Puff_v_air = 20  # m/s
Hyb_v_air = 20  # m/s
MS_v_air = 20  # m/s
wind = 9  # m/s

Puff_t_BZW = (BWZ_dist / (Puff_v_air - wind)) / 60  # min
# print(Puff_t_BZW,'min :time it takes to fly to wind farm for Puffin')
Hyb_t_BZW = (BWZ_dist / (Hyb_v_air - wind)) / 60  # min
# print(Hyb_t_BZW,'min :time it takes to fly to wind farm for Hybrid')
MS_t_BZW = (BWZ_dist / (MS_v_air - wind)) / 60  # min
print(Hyb_t_BZW, 'min :time it takes to fly from shore to BWZ')

# On/In water distance
current = 0.65  # m/s
wind_surface = 5.9  # m/s
nr_unit = 15
SF = 1.6  # Safety_Factor

# On/In water distance Puffin
# TODO: Jan add time for Puffin

# On/In water distance Hybrid
Hyb_t_transfer = 80  # s
Hyb_t_pic = 180  # s per unit
Hyb_t_w = (Hyb_t_pic * nr_unit + Hyb_t_transfer * (nr_unit - 1)) / 60 * SF
print('Hybrid: time spent maneuvering on the water:', Hyb_t_w, 'min')

Hyb_t_total = (Hyb_t_w + Hyb_t_BZW * 2) / 60
print('Hybrid: Total Mission time for', nr_unit, 'units:', Hyb_t_total, 'hours')

# On/In water distance MS
MS_v = 1.1  # m/s
MS_t_1unit = 100 / MS_v
UUV_v = 2.15  # m/s
MS_t_transfer = 80  # s

MS_t_w = (nr_unit * MS_t_1unit + (nr_unit - 1) * MS_t_transfer) / 60 * SF
print('Multi-system: time spent maneuvering on and in the water:', MS_t_w, 'min')

MS_t_total = (MS_t_w + MS_t_BZW * 2) / 60
print('Multi-system: Total Mission time for', nr_unit, 'units :', MS_t_total, 'hours')

# Energy per hop: MS
C_d = 0.038
rho = 1.225
dist_hop = 100
lipo_energy_density = 150  # Wh/kg

MS_W = 21.4 * 9.81  # N
MS_height_hop = 15  # m
MS_v = 5  # m/s
MS_S = 3  # m^2

E_MS = (2.5 * MS_W * MS_height_hop + C_d * 0.5 * MS_v ** 2 * MS_S * dist_hop) / 3600  # Wh
MS_batt_m = E_MS / 150 * nr_unit  # kg
print('Energy for one hop MS:', E_MS, '[Wh]')
print('battery mass for MS for', nr_unit, 'units:', MS_batt_m, 'kg')

# Energy per hop: MS
Hyb_W = 15 * 9.81  # N
Hyb_height_hop = 5  # m
Hyb_v = 5  # m/s
Hyb_S = 1.5  # m^2

E_Hyb = (2.5 * Hyb_W * Hyb_height_hop + C_d * (1 / 2) * Hyb_v ** 2 * Hyb_S * dist_hop) / 3600
Hyb_batt_m = E_MS / 150 * (nr_unit * 2 - 1)  # kg
print('Energy for one hop Hybrid:', E_Hyb, '[Wh]')
print('battery mass for Hybrid for', nr_unit, 'units:', Hyb_batt_m, 'kg')
