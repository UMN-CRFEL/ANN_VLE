import numpy as np

num_species = 4
Species_names = ['N2','C6H14','H2O','C12H26']
SpeciesData = {}

for i in range(num_species):
    SpeciesData[Species_names[i]] = {}

#data for N2
SpeciesData['N2']['W'] = 28.0134 #g/mol
SpeciesData['N2']['Tc'] = 126.1 #K
SpeciesData['N2']['rhoc'] = 11.18 #mol/l
SpeciesData['N2']['Pc'] = 33.978E+5 #Pa
SpeciesData['N2']['omega'] = 0.039
SpeciesData['N2']['thermo'] = {}
SpeciesData['N2']['thermo']['Cp'] = 1005
SpeciesData['N2']['thermo']['Tlow'] = 100
SpeciesData['N2']['thermo']['Tcommon'] = 1000
SpeciesData['N2']['thermo']['Thigh'] = 5000
SpeciesData['N2']['thermo']['lowCpCoeffs'] = np.array([3.298677000E+00,1.408240400E-03,-3.963222000E-06,5.641515000E-09,-2.444854000E-12,-1.020899900E+03,3.950372000E+00])
SpeciesData['N2']['thermo']['highCpCoeffs'] = np.array([2.926640000E+00,1.487976800E-03,-5.684760000E-07,1.009703800E-10,-6.753351000E-15,-9.227977000E+02,5.980528000E+00])


#data for C6H14
SpeciesData['C6H14']['W'] = 86.1754 #g/mol
SpeciesData['C6H14']['Tc'] = 507.6 #K
SpeciesData['C6H14']['rhoc'] = 2.71 #mol/l
SpeciesData['C6H14']['Pc'] = 30.2E+5 #Pa
SpeciesData['C6H14']['omega'] = 0.3003
SpeciesData['C6H14']['thermo'] = {}
SpeciesData['C6H14']['thermo']['Cp'] = 1005
SpeciesData['C6H14']['thermo']['Tlow'] = 100
SpeciesData['C6H14']['thermo']['Tcommon'] = 1200
SpeciesData['C6H14']['thermo']['Thigh'] = 5000
SpeciesData['C6H14']['thermo']['lowCpCoeffs'] = np.array([1.90038413E+00,5.14640294E-02,7.91906421E-06,-3.47484905E-08,1.31381393E-11,-2.29428369E+04,2.05073363E+01])
SpeciesData['C6H14']['thermo']['highCpCoeffs'] = np.array([2.28046470E+01,2.09800119E-02,-3.53074129E-06,-5.46609072E-10,1.47893745E-13,-3.07498496E+04,-9.58469533E+01])

#data for H2O
SpeciesData['H2O']['W'] = 18.0153 #g/mol
SpeciesData['H2O']['Tc'] = 647 #K
SpeciesData['H2O']['rhoc'] = 17.9 #mol/l
SpeciesData['H2O']['Pc'] = 220.64E+5 #Pa
SpeciesData['H2O']['omega'] = 0.344
SpeciesData['H2O']['thermo'] = {}
SpeciesData['H2O']['thermo']['Cp'] = 1005
SpeciesData['H2O']['thermo']['Tlow'] = 100
SpeciesData['H2O']['thermo']['Tcommon'] = 1000
SpeciesData['H2O']['thermo']['Thigh'] = 3500
SpeciesData['H2O']['thermo']['lowCpCoeffs'] = np.array([4.198640560E+00,-2.036434100E-03,6.520402110E-06,-5.487970620E-09,1.771978170E-12,-3.029372670E+04,-8.490322080E-01])
SpeciesData['H2O']['thermo']['highCpCoeffs'] = np.array([3.033992490E+00,2.176918040E-03,-1.640725180E-07,-9.704198700E-11,1.682009920E-14,-3.000429710E+04,4.966770100E+00])

#data for C12H26
SpeciesData['C12H26']['W'] = 170.3348 #g/mol
SpeciesData['C12H26']['Tc'] = 658.2 #K
SpeciesData['C12H26']['rhoc'] = 1.3 #mol/l
SpeciesData['C12H26']['Pc'] = 18.0E+5 #Pa
SpeciesData['C12H26']['omega'] = 0.576385
SpeciesData['C12H26']['thermo'] = {}
SpeciesData['C12H26']['thermo']['Cp'] = 1005
SpeciesData['C12H26']['thermo']['Tlow'] = 100
SpeciesData['C12H26']['thermo']['Tcommon'] = 1590
SpeciesData['C12H26']['thermo']['Thigh'] = 5000
SpeciesData['C12H26']['thermo']['lowCpCoeffs'] = np.array([-2.38265893E+00,1.45739929E-01,-9.14517778E-05,2.85289455E-08,-3.49138416E-12,-4.00920589E+04,4.90709953E+01])
SpeciesData['C12H26']['thermo']['highCpCoeffs'] = np.array([3.85078111E+01,5.63574461E-02,-1.91505499E-05,2.96050890E-09,-1.71263883E-13,-5.48939801E+04,-1.72672880E+02])

SpeciesData['beta'] = np.array([[0,0.02,0.1738,0.19],
                                [0.02,0,0,0],
                                [0.1738,0,0,0],
                                [0.19,0,0,0]])
