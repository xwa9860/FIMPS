#!/usr/bin/python
from pebble import *
class TEST_Reactor:

    def __init__(self):

        import mocup

        self.Pebble = Pebble()
        self.title = ''
        self.T_bulk = 650. + 273.15
        self.Re = 1200.
        self.power = 20.
        # power of test reactor (MWth)

        self.characteristic_radius = 4.27098E+01
        # characteristc radius in (cm)
        self.active_height = 130.3917
        # hieght of active region in (cm)
        self.converging_angle = 45
        # converging angle (deg)
        self.defueling_radius_normalized = 10.
        # radius of defueling channel in terms of pebble radii
        self.wall_effects_normalized = 2
        # thickness of wall effects area in terms of pebble radii
        self.low_porosity = .37
        self.high_porosity = .43
        self.dr_barrel = .3
        # thickness of core barrel (cm)
        self.dr_pressure_vessel = .3
        # thickness of pressure vessel (cm)
        self.buffer_volume_normalized = .2
        # volume of excess volume below active region of the core
        self.defueling_chute_volume_normalized = .1
        # fraction of volume of the defueling chute
        self.lower_plenum_aspect_ratio = 1./4.
        # ratio of diameter of down comer to depth of the lower plenum

        self.shutdown_rods = 12
        # the number of control or shutdown rods - modeled as channels of
        # coolant

        self.mat_steel = mocup.material()
        self.mat_steel.comp['280580'] = 0.0727937298976705
        self.mat_steel.comp['280600'] = 0.0271052041042211
        self.mat_steel.comp['280610'] = 0.00115903506102997
        self.mat_steel.comp['280620'] = 0.00363508684378547
        self.mat_steel.comp['280640'] = 0.000897330790259578
        self.mat_steel.comp['240500'] = 0.000538941136404716
        self.mat_steel.comp['240520'] = 0.00999321476458678
        self.mat_steel.comp['240530'] = 0.00111177014069114
        self.mat_steel.comp['240540'] = 0.00027161824736759
        self.mat_steel.comp['420920'] = 0.00257242129449192
        self.mat_steel.comp['420940'] = 0.00156931418937648
        self.mat_steel.comp['420950'] = 0.00267248625674893
        self.mat_steel.comp['420960'] = 0.00277089989167475
        self.mat_steel.comp['420970'] = 0.00157010125138286
        self.mat_steel.comp['420980'] = 0.00392669589699672
        self.mat_steel.comp['421000'] = 0.00153575631406203
        self.mat_steel.comp['260540'] = 0.000475803634772582
        self.mat_steel.comp['260560'] = 0.00725553665822793
        self.mat_steel.comp['260570'] = 0.000170978438629712
        self.mat_steel.comp['260580'] = 2.13857043928699E-05

        self.mat_graphite = mocup.material()
        self.mat_graphite.comp['60120'] = 1.74/12.
        self.scat = 'grph'

        self.mat_flibe = mocup.material()
        self.mat_flibe.comp['30060'] = 2.44179273755529E-07
        self.mat_flibe.comp['30070'] = 0.0241759215390954
        self.mat_flibe.comp['40090'] = 0.0120880828591846
        self.mat_flibe.comp['90190'] = 0.0483523314367383

    def inputMCNP5(self):

        import math

        # generate geometric parameters for test reactor
        active_volume = math.pi * \
            math.pow(self.characteristic_radius, 2.)*self.active_height
        active_volume += (1./3.)*math.pi*(math.pow(self.characteristic_radius,
                                                   3.)-math.pow((self.defueling_radius_normalized*self.Pebble.r_shell),
                                                                3.))*math.tan(self.converging_angle*math.pi/180.)
        power_density = self.power/(active_volume*1.e-6)

        # generate temperature distribution in test reactor
        self.Pebble.HT(self.T_bulk, self.Re, power_density)

        skeleton = 'Fullcore PB-AHTR\nc ------------------------------------------------------------------------------\nc Cells\nc ------------------------------------------------------------------------------\nc ------------------------------------------------------------------------------\nc TRISO Particles\nc ------------------------------------------------------------------------------\nc\n_seed_title_\nc\nc Depletion Zone 1\n 10  11  -10.5 -1 u=101 tmp=_fuel_temperature_ imp:n=1 $ Kernel 1\n      vol=_seed_kernel_volume_\n11 19 -_seed_matrix_rho_  1 u=101 tmp=8.81648E-08 imp:n=1 $ Coatings and Matrix \n12 0 -11 12 -13 14 -15 16  lat=1            imp:n=1 fill=101 u=81 $ TRISO lattice\nc Depletion Zone 2\n 15  12  -10.5 -1 u=102 tmp=_fuel_temperature_ imp:n=1 $ Kernel 2\n      vol=_seed_kernel_volume_\n16 19 -_seed_matrix_rho_  1 u=102 tmp=8.81648E-08 imp:n=1 $ Coatings and Matrix \n17 0 -11 12 -13 14 -15 16  lat=1            imp:n=1 fill=102 u=82 $ TRISO lattice\nc Depletion Zone 3\n 20  13  -10.5 -1 u=103 tmp=_fuel_temperature_ imp:n=1 $ Kernel 3\n      vol=_seed_kernel_volume_\n21 19 -_seed_matrix_rho_  1 u=103 tmp=8.81648E-08 imp:n=1 $ Coatings and Matrix \n22 0 -11 12 -13 14 -15 16  lat=1            imp:n=1 fill=103 u=83 $ TRISO lattice\nc Depletion Zone 4\n 25  14  -10.5 -1 u=104 tmp=_fuel_temperature_ imp:n=1 $ Kernel 4\n      vol=_seed_kernel_volume_\n26 19 -_seed_matrix_rho_  1 u=104 tmp=8.81648E-08 imp:n=1 $ Coatings and Matrix \n27 0 -11 12 -13 14 -15 16  lat=1            imp:n=1 fill=104 u=84 $ TRISO \nc\nc ------------------------------------------------------------------------------\nc Pebbles\nc ------------------------------------------------------------------------------\nc\nc Seed Pebbles \nc\n110 1 -_seed_core_rho_      -110             u=11 imp:n=1 $ Depletion Region 1\n          tmp=_seed_core_temperature_\n111 0                     -111 110 fill=81 u=11 imp:n=1\n112 1 -_seed_shell_rho_     -112 111         u=11 imp:n=1\n          tmp=_seed_shell_temperature_\n113 1 -_seed_core_rho_      -113             u=11 imp:n=1 $ Depletion Region 1\n          tmp=_seed_core_temperature_\n114 0                     -114 113 fill=81 u=11 imp:n=1\n115 1 -_seed_shell_rho_     -115 114         u=11 imp:n=1\n          tmp=_seed_shell_temperature_\n116 1 -_seed_core_rho_      -116             u=11 imp:n=1 $ Depletion Region 1\n          tmp=_seed_core_temperature_\n117 0                     -117 116 fill=81 u=11 imp:n=1\n118 1 -_seed_shell_rho_     -118 117         u=11 imp:n=1\n          tmp=_seed_shell_temperature_\n119 1 -_seed_core_rho_      -119             u=11 imp:n=1 $ Depletion Region 1\n          tmp=_seed_core_temperature_\n120 0                     -120 119 fill=81 u=11 imp:n=1\n121 1 -_seed_shell_rho_     -121 120         u=11 imp:n=1\n          tmp=_seed_shell_temperature_\nc\n123 1 -_seed_core_rho_      -130             u=11 imp:n=1 $ Depletion Region 1\n          tmp=_seed_core_temperature_\n124 0                     -131 130 fill=81 u=11 imp:n=1\n125 1 -_seed_shell_rho_     -132 131         u=11 imp:n=1\n          tmp=_seed_shell_temperature_\n126 1 -_seed_core_rho_      -133             u=11 imp:n=1 $ Depletion Region 1\n          tmp=_seed_core_temperature_\n127 0                     -134 133 fill=81 u=11 imp:n=1\n128 1 -_seed_shell_rho_     -135 134         u=11 imp:n=1\n          tmp=_seed_shell_temperature_\n129 1 -_seed_core_rho_      -136             u=11 imp:n=1 $ Depletion Region 1\n          tmp=_seed_core_temperature_\n130 0                     -137 136 fill=81 u=11 imp:n=1\n131 1 -_seed_shell_rho_     -138 137         u=11 imp:n=1\n          tmp=_seed_shell_temperature_\n132 1 -_seed_core_rho_      -139             u=11 imp:n=1 $ Depletion Region 1\n          tmp=_seed_core_temperature_\n133 0                     -140 139 fill=81 u=11 imp:n=1\n134 1 -_seed_shell_rho_     -141 140         u=11 imp:n=1\n          tmp=_seed_shell_temperature_\nc\n135 1 -_seed_core_rho_      -150             u=11 imp:n=1 $ Depletion Region 2\n          tmp=_seed_core_temperature_\n136 0                     -151 150 fill=82 u=11 imp:n=1\n137 1 -_seed_shell_rho_     -152 151         u=11 imp:n=1\n          tmp=_seed_shell_temperature_\nc\n138 1 -_seed_core_rho_      -153             u=11 imp:n=1 $ Depletion Region 2\n          tmp=_seed_core_temperature_\n139 0                     -154 153 fill=82 u=11 imp:n=1\n140 1 -_seed_shell_rho_     -155 154         u=11 imp:n=1\n          tmp=_seed_shell_temperature_\nc\n141 1 -_seed_core_rho_      -160             u=11 imp:n=1 $ Depletion Region 3\n          tmp=_seed_core_temperature_\n142 0                     -161 160 fill=83 u=11 imp:n=1\n143 1 -_seed_shell_rho_     -162 161         u=11 imp:n=1\n          tmp=_seed_shell_temperature_\nc\n144 1 -_seed_core_rho_      -163             u=11 imp:n=1 $ Depletion Region 3\n          tmp=_seed_core_temperature_\n145 0                     -164 163 fill=83 u=11 imp:n=1\n146 1 -_seed_shell_rho_     -165 164         u=11 imp:n=1\n          tmp=_seed_shell_temperature_\nc\n147 1 -_seed_core_rho_      -170             u=11 imp:n=1 $ Depletion Region 4\n          tmp=_seed_core_temperature_\n148 0                     -171 170 fill=84 u=11 imp:n=1\n149 1 -_seed_shell_rho_     -172 171         u=11 imp:n=1\n          tmp=_seed_shell_temperature_\nc\n150 1 -_seed_core_rho_      -173             u=11 imp:n=1 $ Depletion Region 4\n          tmp=_seed_core_temperature_\n151 0                     -174 173 fill=84 u=11 imp:n=1\n152 1 -_seed_shell_rho_     -175 174         u=11 imp:n=1\n          tmp=_seed_shell_temperature_\nc\n153 30 -_coolant_density_ 112 115 118 121 \n                   132 135 138 141\n                   152 155 162 165 172 175 u=11 imp:n=1 $ Coolant\n          tmp=_coolant_temperature_\n154 0              -100 101 -102 103 -104 105 fill=11 imp:n=1 u=1 lat=1 $ FCC Unit Cell\nc\nc Seed Pebbles (converging)\nc\n210 1 _seed_core_rho_      -210             u=21 imp:n=1 $ Depletion Region 1\n          tmp=_seed_core_temperature_\n211 0                     -211 210 fill=81 u=21 imp:n=1\n212 1 _seed_shell_rho_     -212 211         u=21 imp:n=1\n          tmp=_seed_shell_temperature_\n213 1 _seed_core_rho_      -213             u=21 imp:n=1 $ Depletion Region 1\n          tmp=_seed_core_temperature_\n214 0                     -214 213 fill=81 u=21 imp:n=1\n215 1 _seed_shell_rho_     -215 214         u=21 imp:n=1\n          tmp=_seed_shell_temperature_\n216 1 _seed_core_rho_      -216             u=21 imp:n=1 $ Depletion Region 1\n          tmp=_seed_core_temperature_\n217 0                     -217 216 fill=81 u=21 imp:n=1\n218 1 _seed_shell_rho_     -218 217         u=21 imp:n=1\n          tmp=_seed_shell_temperature_\n219 1 _seed_core_rho_      -219             u=21 imp:n=1 $ Depletion Region 1\n          tmp=_seed_core_temperature_\n220 0                     -220 219 fill=81 u=21 imp:n=1\n221 1 _seed_shell_rho_     -221 220         u=21 imp:n=1\n          tmp=_seed_shell_temperature_\nc\n223 1 _seed_core_rho_      -230             u=21 imp:n=1 $ Depletion Region 1\n          tmp=_seed_core_temperature_\n224 0                     -231 230 fill=81 u=21 imp:n=1\n225 1 _seed_shell_rho_     -232 231         u=21 imp:n=1\n          tmp=_seed_shell_temperature_\n226 1 _seed_core_rho_      -233             u=21 imp:n=1 $ Depletion Region 1\n          tmp=_seed_core_temperature_\n227 0                     -234 233 fill=81 u=21 imp:n=1\n228 1 _seed_shell_rho_     -235 234         u=21 imp:n=1\n          tmp=_seed_shell_temperature_\n229 1 _seed_core_rho_      -236             u=21 imp:n=1 $ Depletion Region 1\n          tmp=_seed_core_temperature_\n230 0                     -237 236 fill=81 u=21 imp:n=1\n231 1 _seed_shell_rho_     -238 237         u=21 imp:n=1\n          tmp=_seed_shell_temperature_\n232 1 _seed_core_rho_      -239             u=21 imp:n=1 $ Depletion Region 1\n          tmp=_seed_core_temperature_\n233 0                     -240 239 fill=81 u=21 imp:n=1\n234 1 _seed_shell_rho_     -241 240         u=21 imp:n=1\n          tmp=_seed_shell_temperature_\nc\n235 1 _seed_core_rho_      -250             u=21 imp:n=1 $ Depletion Region 2\n          tmp=_seed_core_temperature_\n236 0                     -251 250 fill=82 u=21 imp:n=1\n237 1 _seed_shell_rho_     -252 251         u=21 imp:n=1\n          tmp=_seed_shell_temperature_\nc\n238 1 _seed_core_rho_      -253             u=21 imp:n=1 $ Depletion Region 2\n          tmp=_seed_core_temperature_\n239 0                     -254 253 fill=82 u=21 imp:n=1\n240 1 _seed_shell_rho_     -255 254         u=21 imp:n=1\n          tmp=_seed_shell_temperature_\nc\n241 1 _seed_core_rho_      -260             u=21 imp:n=1 $ Depletion Region 3\n          tmp=_seed_core_temperature_\n242 0                     -261 260 fill=83 u=21 imp:n=1\n243 1 _seed_shell_rho_     -262 261         u=21 imp:n=1\n          tmp=_seed_shell_temperature_\nc\n244 1 _seed_core_rho_      -263             u=21 imp:n=1 $ Depletion Region 3\n          tmp=_seed_core_temperature_\n245 0                     -264 263 fill=83 u=21 imp:n=1\n246 1 _seed_shell_rho_     -265 264         u=21 imp:n=1\n          tmp=_seed_shell_temperature_\nc\n247 1 _seed_core_rho_      -270             u=21 imp:n=1 $ Depletion Region 4\n          tmp=_seed_core_temperature_\n248 0                     -271 270 fill=84 u=21 imp:n=1\n249 1 _seed_shell_rho_     -272 271         u=21 imp:n=1\n          tmp=_seed_shell_temperature_\nc\n250 1 _seed_core_rho_      -273             u=21 imp:n=1 $ Depletion Region 4\n          tmp=_seed_core_temperature_\n251 0                     -274 273 fill=84 u=21 imp:n=1\n252 1 _seed_shell_rho_     -275 274         u=21 imp:n=1\n          tmp=_seed_shell_temperature_\nc\n253 30 -_coolant_density_ 212 215 218 221 \n                   232 235 238 241\n                   252 255 262 265 272 275 u=21 imp:n=1 $ Coolant\n          tmp=_coolant_temperature_\n254 0              -200 201 -202 203 -204 205 fill=21 imp:n=1 u=2 lat=1 $ FCC Unit Cell\nc\nc Inner Graphite Reflector\nc\n_Outer_Reflector_\n_Shutdown_Channels_\n901 31 -8.86       -905  904 -805 812 imp:n=1        $ core barrel\n          tmp=7.52393E-08\n902 30 -1.9872 -906  905 -805 812 imp:n=1        $ coolant in downcomer \n          tmp=7.52393E-08\n903 31 -8.86       -907  906 -805 812 imp:n=1        $ reactor pressure vessel\n          tmp=7.52393E-08\nc\nc Active Region\nc\n920 30 -1.9872 -923      -802 812 imp:n=1        $ coolant below pebble free surface\n          tmp=7.52393E-08\n921 13 -10.5 -922      -803 802 imp:n=1 fill=1 $ active region \n922 13 -10.5 -923 922  -803 802 imp:n=1 fill=2 $ active region (wall effects)\nc\nc Contraction Region\nc\n931 13 -10.5 -933      -804 803 imp:n=1 fill=2 $ pebbles in converging region\n934 1  -1.74 -923  933 -804 803 imp:n=1        $ outer radial graphite reflector lower level\n          tmp=8.38563E-08\nc\nc Upper Level\nc\n940 13 -10.5 -943     -805  804 imp:n=1 fill=2 $ pebbles in defueling chute\n950 1  -1.74 -923  943 -805 804 imp:n=1        $ outer radial graphite reflector lower level\n          tmp=8.38563E-08\nc\nc Lower Plenum\nc\n910 30 -1.9872         -812 -801      imp:n=1   $ Lower Coolant Plenum\n          tmp=7.52393E-08\n911 31 -8.86           -812  801 -800 imp:n=1   $ Lower Reactor Pressure Vessel\n          tmp=7.52393E-08\nc Vacuum Boundary Conditions\nc\n999 0              907:805:(800 -812) imp:n=0        $ void\n\nc ------------------------------------------------------------------------------\nc Surfaces\nc ------------------------------------------------------------------------------\nc TRISO surfaces\nc Seed\n1   so    _seed_kernel_r_    $  Kernel\nc Blanket\nc   Seed Fuel Particles Lattice\n11  px    _seed_TRISO_hp_\n12  px   -_seed_TRISO_hp_\n13  py    _seed_TRISO_hp_\n14  py   -_seed_TRISO_hp_\n15  pz    _seed_TRISO_hp_\n16  pz   -_seed_TRISO_hp_\nc\nc Seed Pebble Surfaces\n110 s  _seed_pebble_hp1_  _seed_pebble_hp1_  _seed_pebble_hp1_ _seed_core_r_ \n111 s  _seed_pebble_hp1_  _seed_pebble_hp1_  _seed_pebble_hp1_ _seed_active_r_\n112 s  _seed_pebble_hp1_  _seed_pebble_hp1_  _seed_pebble_hp1_ _seed_shell_r_\nc\n113 s -_seed_pebble_hp1_  _seed_pebble_hp1_  _seed_pebble_hp1_ _seed_core_r_ \n114 s -_seed_pebble_hp1_  _seed_pebble_hp1_  _seed_pebble_hp1_ _seed_active_r_\n115 s -_seed_pebble_hp1_  _seed_pebble_hp1_  _seed_pebble_hp1_ _seed_shell_r_\nc\n116 s -_seed_pebble_hp1_ -_seed_pebble_hp1_  _seed_pebble_hp1_ _seed_core_r_ \n117 s -_seed_pebble_hp1_ -_seed_pebble_hp1_  _seed_pebble_hp1_ _seed_active_r_\n118 s -_seed_pebble_hp1_ -_seed_pebble_hp1_  _seed_pebble_hp1_ _seed_shell_r_\nc\n119 s  _seed_pebble_hp1_ -_seed_pebble_hp1_  _seed_pebble_hp1_ _seed_core_r_ \n120 s  _seed_pebble_hp1_ -_seed_pebble_hp1_  _seed_pebble_hp1_ _seed_active_r_\n121 s  _seed_pebble_hp1_ -_seed_pebble_hp1_  _seed_pebble_hp1_ _seed_shell_r_\nc\n130 s  _seed_pebble_hp1_  _seed_pebble_hp1_ -_seed_pebble_hp1_ _seed_core_r_ \n131 s  _seed_pebble_hp1_  _seed_pebble_hp1_ -_seed_pebble_hp1_ _seed_active_r_ \n132 s  _seed_pebble_hp1_  _seed_pebble_hp1_ -_seed_pebble_hp1_ _seed_shell_r_\nc\n133 s -_seed_pebble_hp1_  _seed_pebble_hp1_ -_seed_pebble_hp1_ _seed_core_r_\n134 s -_seed_pebble_hp1_  _seed_pebble_hp1_ -_seed_pebble_hp1_ _seed_active_r_\n135 s -_seed_pebble_hp1_  _seed_pebble_hp1_ -_seed_pebble_hp1_ _seed_shell_r_\nc\n136 s -_seed_pebble_hp1_ -_seed_pebble_hp1_ -_seed_pebble_hp1_ _seed_core_r_\n137 s -_seed_pebble_hp1_ -_seed_pebble_hp1_ -_seed_pebble_hp1_ _seed_active_r_\n138 s -_seed_pebble_hp1_ -_seed_pebble_hp1_ -_seed_pebble_hp1_ _seed_shell_r_\nc\n139 s  _seed_pebble_hp1_ -_seed_pebble_hp1_ -_seed_pebble_hp1_ _seed_core_r_\n140 s  _seed_pebble_hp1_ -_seed_pebble_hp1_ -_seed_pebble_hp1_ _seed_active_r_\n141 s  _seed_pebble_hp1_ -_seed_pebble_hp1_ -_seed_pebble_hp1_ _seed_shell_r_\nc\n150 s  _seed_pebble_hp1_  0           0          _seed_core_r_\n151 s  _seed_pebble_hp1_  0           0          _seed_active_r_\n152 s  _seed_pebble_hp1_  0           0          _seed_shell_r_\nc\n153 s -_seed_pebble_hp1_  0           0          _seed_core_r_\n154 s -_seed_pebble_hp1_  0           0          _seed_active_r_\n155 s -_seed_pebble_hp1_  0           0          _seed_shell_r_\nc\n160 s   0          _seed_pebble_hp1_  0          _seed_core_r_\n161 s   0          _seed_pebble_hp1_  0          _seed_active_r_\n162 s   0          _seed_pebble_hp1_  0          _seed_shell_r_\nc\n163 s   0         -_seed_pebble_hp1_  0          _seed_core_r_\n164 s   0         -_seed_pebble_hp1_  0          _seed_active_r_\n165 s   0         -_seed_pebble_hp1_  0          _seed_shell_r_\nc\n170 s   0           0          _seed_pebble_hp1_  _seed_core_r_\n171 s   0           0          _seed_pebble_hp1_  _seed_active_r_\n172 s   0           0          _seed_pebble_hp1_  _seed_shell_r_\nc\n173 s   0           0         -_seed_pebble_hp1_  _seed_core_r_\n174 s   0           0         -_seed_pebble_hp1_  _seed_active_r_\n175 s   0           0         -_seed_pebble_hp1_  _seed_shell_r_\nc\nc Seed Pebble Surfaces\n210 s  _seed_pebble_hp2_  _seed_pebble_hp2_  _seed_pebble_hp2_ _seed_core_r_ \n211 s  _seed_pebble_hp2_  _seed_pebble_hp2_  _seed_pebble_hp2_ _seed_active_r_\n212 s  _seed_pebble_hp2_  _seed_pebble_hp2_  _seed_pebble_hp2_ _seed_shell_r_\nc\n213 s -_seed_pebble_hp2_  _seed_pebble_hp2_  _seed_pebble_hp2_ _seed_core_r_ \n214 s -_seed_pebble_hp2_  _seed_pebble_hp2_  _seed_pebble_hp2_ _seed_active_r_\n215 s -_seed_pebble_hp2_  _seed_pebble_hp2_  _seed_pebble_hp2_ _seed_shell_r_\nc\n216 s -_seed_pebble_hp2_ -_seed_pebble_hp2_  _seed_pebble_hp2_ _seed_core_r_ \n217 s -_seed_pebble_hp2_ -_seed_pebble_hp2_  _seed_pebble_hp2_ _seed_active_r_\n218 s -_seed_pebble_hp2_ -_seed_pebble_hp2_  _seed_pebble_hp2_ _seed_shell_r_\nc\n219 s  _seed_pebble_hp2_ -_seed_pebble_hp2_  _seed_pebble_hp2_ _seed_core_r_ \n220 s  _seed_pebble_hp2_ -_seed_pebble_hp2_  _seed_pebble_hp2_ _seed_active_r_\n221 s  _seed_pebble_hp2_ -_seed_pebble_hp2_  _seed_pebble_hp2_ _seed_shell_r_\nc\n230 s  _seed_pebble_hp2_  _seed_pebble_hp2_ -_seed_pebble_hp2_ _seed_core_r_ \n231 s  _seed_pebble_hp2_  _seed_pebble_hp2_ -_seed_pebble_hp2_ _seed_active_r_ \n232 s  _seed_pebble_hp2_  _seed_pebble_hp2_ -_seed_pebble_hp2_ _seed_shell_r_\nc\n233 s -_seed_pebble_hp2_  _seed_pebble_hp2_ -_seed_pebble_hp2_ _seed_core_r_\n234 s -_seed_pebble_hp2_  _seed_pebble_hp2_ -_seed_pebble_hp2_ _seed_active_r_\n235 s -_seed_pebble_hp2_  _seed_pebble_hp2_ -_seed_pebble_hp2_ _seed_shell_r_\nc\n236 s -_seed_pebble_hp2_ -_seed_pebble_hp2_ -_seed_pebble_hp2_ _seed_core_r_\n237 s -_seed_pebble_hp2_ -_seed_pebble_hp2_ -_seed_pebble_hp2_ _seed_active_r_\n238 s -_seed_pebble_hp2_ -_seed_pebble_hp2_ -_seed_pebble_hp2_ _seed_shell_r_\nc\n239 s  _seed_pebble_hp2_ -_seed_pebble_hp2_ -_seed_pebble_hp2_ _seed_core_r_\n240 s  _seed_pebble_hp2_ -_seed_pebble_hp2_ -_seed_pebble_hp2_ _seed_active_r_\n241 s  _seed_pebble_hp2_ -_seed_pebble_hp2_ -_seed_pebble_hp2_ _seed_shell_r_\nc\n250 s  _seed_pebble_hp2_  0           0          _seed_core_r_\n251 s  _seed_pebble_hp2_  0           0          _seed_active_r_\n252 s  _seed_pebble_hp2_  0           0          _seed_shell_r_\nc\n253 s -_seed_pebble_hp2_  0           0          _seed_core_r_\n254 s -_seed_pebble_hp2_  0           0          _seed_active_r_\n255 s -_seed_pebble_hp2_  0           0          _seed_shell_r_\nc\n260 s   0          _seed_pebble_hp2_  0          _seed_core_r_\n261 s   0          _seed_pebble_hp2_  0          _seed_active_r_\n262 s   0          _seed_pebble_hp2_  0          _seed_shell_r_\nc\n263 s   0         -_seed_pebble_hp2_  0          _seed_core_r_\n264 s   0         -_seed_pebble_hp2_  0          _seed_active_r_\n265 s   0         -_seed_pebble_hp2_  0          _seed_shell_r_\nc\n270 s   0           0          _seed_pebble_hp2_  _seed_core_r_\n271 s   0           0          _seed_pebble_hp2_  _seed_active_r_\n272 s   0           0          _seed_pebble_hp2_  _seed_shell_r_\nc\n273 s   0           0         -_seed_pebble_hp2_  _seed_core_r_\n274 s   0           0         -_seed_pebble_hp2_  _seed_active_r_\n275 s   0           0         -_seed_pebble_hp2_  _seed_shell_r_\nc\nc Seed Unit Cell Surfaces\n100 px  _seed_pebble_hp1_\n101 px -_seed_pebble_hp1_\n102 py  _seed_pebble_hp1_\n103 py -_seed_pebble_hp1_\n104 pz  _seed_pebble_hp1_\n105 pz -_seed_pebble_hp1_\nc\nc Seed Unit Cell Surfaces\n200 px  _seed_pebble_hp2_\n201 px -_seed_pebble_hp2_\n202 py  _seed_pebble_hp2_\n203 py -_seed_pebble_hp2_\n204 pz  _seed_pebble_hp2_\n205 pz -_seed_pebble_hp2_\nc\n_SR_Channel_Surfaces_\nc\n904 cz _OGRR_\n905 cz _core_barrel_r_\n906 cz _down_comer_r_\n907 cz _pressure_vessel_r_\nc\n922 cz _boundary_layer_R_\n923 cz _active_R_\nc\n933 kz _contraction_b_ _contraction_tan_ -1\nc \n943 cz _exit_R_\nc\nc axial subdivision\nc\n800 sq _plenum1_ 0 0 0 -1 0 0 0 $ Lower Surface of Reactor Pressure Vessel\n801 sq _plenum2_ 0 0 0 -1 0 0 0 $ Lower Surface of Coolant Plenum\n812 pz 0\n802 pz _free_surface_\n803 pz _active_Z_\n804 pz _contraction_Z_\n805 pz _exit_Z_\n\nc ------------------------------------------------------------------------------\nc Data\nc ------------------------------------------------------------------------------\nc   Fuel Mix (800 C - Scattering Kernel at 727 C)\nm11        6000._fuel_XS_ 1.182344e-02\n           8016._fuel_XS_ 3.547031e-02\n          90232._fuel_XS_ 2.239743E-13\n          90233._fuel_XS_ 1.522581E-18\n          91233._fuel_XS_ 1.572958E-14\n          92233._fuel_XS_ 5.571045E-12\n          92234._fuel_XS_ 6.436331E-09\n          92235._fuel_XS_ 4.086884E-03\n          92236._fuel_XS_ 1.279335E-04\n          92237._fuel_XS_ 8.067304E-07\n          92238._fuel_XS_ 1.869822E-02\n          92239._fuel_XS_ 1.459012E-07\n          93236._fuel_XS_ 5.939347E-14\n          93237._fuel_XS_ 1.525899E-06\n          93238._fuel_XS_ 1.351391E-08\n          93239._fuel_XS_ 1.964019E-05\n          94236._fuel_XS_ 9.276975E-14\n          94237._fuel_XS_ 6.866350E-14\n          94238._fuel_XS_ 6.056534E-08\n          94239._fuel_XS_ 1.220318E-04\n          94240._fuel_XS_ 1.481693E-05\n          94241._fuel_XS_ 4.694728E-06\n          94242._fuel_XS_ 2.699156E-07\n          94243._fuel_XS_ 2.583281E-10\n          94244._fuel_XS_ 3.699372E-12\n          95241._fuel_XS_ 7.441561E-09\n          95642._fuel_XS_ 1.758807E-10\n          95242._fuel_XS_ 6.405703E-11\n          95243._fuel_XS_ 9.252288E-09\n          95244._fuel_XS_ 1.393436E-12\n          96242._fuel_XS_ 5.926078E-10\n          96243._fuel_XS_ 2.274614E-12\n          96244._fuel_XS_ 4.561495E-10\n          96245._fuel_XS_ 6.826465E-12\n          96246._fuel_XS_ 8.451498E-14\n          96247._fuel_XS_ 1.774159E-16\n          96248._fuel_XS_ 1.941569E-18\n          97249._fuel_XS_ 6.556835E-21\n          98249._fuel_XS_ 6.628736E-23\n          98250._fuel_XS_ 1.005694E-21\n          35081._fuel_XS_ 1.198639E-06\n          36083._fuel_XS_ 2.869652E-06\n          36084._fuel_XS_ 5.816297E-06\n          36086._fuel_XS_ 1.107374E-05\n          37085._fuel_XS_ 5.505316E-06\n          37087._fuel_XS_ 1.413340E-05\n          38088._fuel_XS_ 2.045563E-05\n          38089._fuel_XS_ 1.998812E-05\n          38090._fuel_XS_ 3.142291E-05\n          39089._fuel_XS_ 6.309500E-06\n          39091._fuel_XS_ 2.534601E-05\n          40091._fuel_XS_ 6.768296E-06\n          40092._fuel_XS_ 4.282915E-05\n          40093._fuel_XS_ 3.533351E-05\n          40094._fuel_XS_ 3.505424E-05\n          40095._fuel_XS_ 2.896422E-05\n          40096._fuel_XS_ 3.552715E-05\n          42095._fuel_XS_ 1.796686E-06\n          42097._fuel_XS_ 3.315333E-05\n          42098._fuel_XS_ 3.344649E-05\n          42099._fuel_XS_ 3.968617E-06\n          42100._fuel_XS_ 3.618368E-05\n          43099._fuel_XS_ 2.936539E-05\n          44101._fuel_XS_ 2.912392E-05\n          44102._fuel_XS_ 2.513771E-05\n          44103._fuel_XS_ 1.377158E-05\n          44104._fuel_XS_ 1.227338E-05\n          44105._fuel_XS_ 6.705035E-08\n          44106._fuel_XS_ 3.610884E-06\n          45103._fuel_XS_ 5.177439E-06\n          45105._fuel_XS_ 4.058725E-07\n          46104._fuel_XS_ 3.231088E-07\n          46105._fuel_XS_ 5.367453E-06\n          46106._fuel_XS_ 1.754949E-06\n          46107._fuel_XS_ 2.147861E-06\n          46108._fuel_XS_ 1.274244E-06\n          46110._fuel_XS_ 4.024549E-07\n          47109._fuel_XS_ 6.889880E-07\n          48110._fuel_XS_ 4.860286E-08\n          48111._fuel_XS_ 1.718458E-07\n          48113._fuel_XS_ 1.165852E-08\n          48114._fuel_XS_ 2.399130E-07\n          49115._fuel_XS_ 7.316891E-08\n          52130._fuel_XS_ 1.000000E-30\n          53127._fuel_XS_ 7.455139E-07\n          53129._fuel_XS_ 4.076083E-06\n          54131._fuel_XS_ 1.102128E-05\n          54132._fuel_XS_ 2.235732E-05\n          54134._fuel_XS_ 4.399099E-05\n          54135._fuel_XS_ 6.329559E-08\n          54136._fuel_XS_ 7.148401E-05\n          55133._fuel_XS_ 2.827298E-05\n          55134._fuel_XS_ 6.614927E-07\n          55135._fuel_XS_ 3.763328E-06\n          55137._fuel_XS_ 3.501798E-05\n          56138._fuel_XS_ 3.884912E-05\n          56140._fuel_XS_ 1.459397E-05\n          57139._fuel_XS_ 3.648224E-05\n          58141._fuel_XS_ 2.202250E-05\n          58142._fuel_XS_ 3.351438E-05\n          58143._fuel_XS_ 1.981300E-06\n          59141._fuel_XS_ 1.099813E-05\n          59143._fuel_XS_ 1.355711E-05\n          60143._fuel_XS_ 1.670550E-05\n          60144._fuel_XS_ 2.522952E-06\n          60145._fuel_XS_ 2.136752E-05\n          60146._fuel_XS_ 1.759809E-05\n          60147._fuel_XS_ 4.659857E-06\n          60148._fuel_XS_ 1.016186E-05\n          60150._fuel_XS_ 3.836309E-06\n          61147._fuel_XS_ 6.886872E-06\n          61148._fuel_XS_ 1.556989E-07\n          61548._fuel_XS_ 4.434201E-08\n          61149._fuel_XS_ 6.082764E-07\n          62147._fuel_XS_ 9.020074E-08\n          62149._fuel_XS_ 3.477574E-07\n          62150._fuel_XS_ 5.469674E-06\n          62151._fuel_XS_ 1.030612E-06\n          62152._fuel_XS_ 2.692367E-06\n          62153._fuel_XS_ 1.485319E-07\n          62154._fuel_XS_ 5.646419E-07\n          63153._fuel_XS_ 1.211060E-06\n          63154._fuel_XS_ 9.311692E-08\n          63155._fuel_XS_ 1.004228E-07\n          63156._fuel_XS_ 1.514866E-07\n          64155._fuel_XS_ 1.908087E-10\n          64156._fuel_XS_ 1.448828E-07\n          64157._fuel_XS_ 4.434819E-09\n          64158._fuel_XS_ 1.170789E-07\nc  \nm12        6000._fuel_XS_ 1.182344e-02\n           8016._fuel_XS_ 3.547031e-02\n          90232._fuel_XS_ 1.394671E-12\n          90233._fuel_XS_ 9.518447E-18\n          91233._fuel_XS_ 1.906775E-13\n          92233._fuel_XS_ 1.210829E-11\n          92234._fuel_XS_ 1.622178E-08\n          92235._fuel_XS_ 2.992394E-03\n          92236._fuel_XS_ 3.290182E-04\n          92237._fuel_XS_ 2.440635E-06\n          92238._fuel_XS_ 1.831094E-02\n          92239._fuel_XS_ 1.436948E-07\n          93236._fuel_XS_ 7.108439E-13\n          93237._fuel_XS_ 1.093487E-05\n          93238._fuel_XS_ 1.056148E-07\n          93239._fuel_XS_ 2.055978E-05\n          94236._fuel_XS_ 1.432859E-12\n          94237._fuel_XS_ 4.628999E-13\n          94238._fuel_XS_ 1.048356E-06\n          94239._fuel_XS_ 2.559519E-04\n          94240._fuel_XS_ 5.936493E-05\n          94241._fuel_XS_ 3.616747E-05\n          94242._fuel_XS_ 4.791085E-06\n          94243._fuel_XS_ 4.690254E-09\n          94244._fuel_XS_ 1.581908E-10\n          95241._fuel_XS_ 1.183441E-07\n          95642._fuel_XS_ 5.184383E-09\n          95242._fuel_XS_ 1.054759E-09\n          95243._fuel_XS_ 3.568993E-07\n          95244._fuel_XS_ 5.489578E-11\n          96242._fuel_XS_ 2.137292E-08\n          96243._fuel_XS_ 1.681119E-10\n          96244._fuel_XS_ 3.721591E-08\n          96245._fuel_XS_ 1.062860E-09\n          96246._fuel_XS_ 2.781395E-11\n          96247._fuel_XS_ 1.161531E-13\n          96248._fuel_XS_ 2.602645E-15\n          97249._fuel_XS_ 1.513246E-17\n          98249._fuel_XS_ 2.861011E-19\n          98250._fuel_XS_ 4.140578E-18\n          35081._fuel_XS_ 3.324822E-06\n          36083._fuel_XS_ 7.559134E-06\n          36084._fuel_XS_ 1.632053E-05\n          36086._fuel_XS_ 3.013223E-05\n          37085._fuel_XS_ 1.512706E-05\n          37087._fuel_XS_ 3.847264E-05\n          38088._fuel_XS_ 5.553302E-05\n          38089._fuel_XS_ 3.830754E-05\n          38090._fuel_XS_ 8.517845E-05\n          39089._fuel_XS_ 3.305612E-05\n          39091._fuel_XS_ 5.160390E-05\n          40091._fuel_XS_ 3.708090E-05\n          40092._fuel_XS_ 1.185987E-04\n          40093._fuel_XS_ 9.874096E-05\n          40094._fuel_XS_ 9.763775E-05\n          40095._fuel_XS_ 6.126198E-05\n          40096._fuel_XS_ 9.958187E-05\n          42095._fuel_XS_ 1.774545E-05\n          42097._fuel_XS_ 9.576307E-05\n          42098._fuel_XS_ 9.505332E-05\n          42099._fuel_XS_ 3.792566E-06\n          42100._fuel_XS_ 1.033081E-04\n          43099._fuel_XS_ 8.883524E-05\n          44101._fuel_XS_ 8.290260E-05\n          44102._fuel_XS_ 7.424203E-05\n          44103._fuel_XS_ 2.847202E-05\n          44104._fuel_XS_ 3.949793E-05\n          44105._fuel_XS_ 8.541761E-08\n          44106._fuel_XS_ 1.386107E-05\n          45103._fuel_XS_ 2.636898E-05\n          45105._fuel_XS_ 5.301878E-07\n          46104._fuel_XS_ 3.765796E-06\n          46105._fuel_XS_ 1.923902E-05\n          46106._fuel_XS_ 7.208345E-06\n          46107._fuel_XS_ 9.466758E-06\n          46108._fuel_XS_ 6.069032E-06\n          46110._fuel_XS_ 1.845983E-06\n          47109._fuel_XS_ 3.309932E-06\n          48110._fuel_XS_ 5.132540E-07\n          48111._fuel_XS_ 8.673683E-07\n          48113._fuel_XS_ 1.566864E-08\n          48114._fuel_XS_ 8.247829E-07\n          49115._fuel_XS_ 2.012776E-07\n          52130._fuel_XS_ 1.000000E-30\n          53127._fuel_XS_ 2.693293E-06\n          53129._fuel_XS_ 1.270155E-05\n          54131._fuel_XS_ 3.747744E-05\n          54132._fuel_XS_ 7.431147E-05\n          54134._fuel_XS_ 1.250791E-04\n          54135._fuel_XS_ 5.621731E-08\n          54136._fuel_XS_ 2.051967E-04\n          55133._fuel_XS_ 9.367238E-05\n          55134._fuel_XS_ 5.187854E-06\n          55135._fuel_XS_ 1.094722E-05\n          55137._fuel_XS_ 9.974388E-05\n          56138._fuel_XS_ 1.090941E-04\n          56140._fuel_XS_ 1.816436E-05\n          57139._fuel_XS_ 1.027218E-04\n          58141._fuel_XS_ 3.849038E-05\n          58142._fuel_XS_ 9.507646E-05\n          58143._fuel_XS_ 1.748777E-06\n          59141._fuel_XS_ 5.429249E-05\n          59143._fuel_XS_ 1.765595E-05\n          60143._fuel_XS_ 6.566093E-05\n          60144._fuel_XS_ 1.753406E-05\n          60145._fuel_XS_ 5.867832E-05\n          60146._fuel_XS_ 5.091266E-05\n          60147._fuel_XS_ 5.461265E-06\n          60148._fuel_XS_ 2.930830E-05\n          60150._fuel_XS_ 1.126892E-05\n          61147._fuel_XS_ 2.191141E-05\n          61148._fuel_XS_ 5.833424E-07\n          61548._fuel_XS_ 1.530990E-07\n          61149._fuel_XS_ 7.063616E-07\n          62147._fuel_XS_ 7.009690E-07\n          62149._fuel_XS_ 4.171437E-07\n          62150._fuel_XS_ 1.839040E-05\n          62151._fuel_XS_ 1.672633E-06\n          62152._fuel_XS_ 8.309547E-06\n          62153._fuel_XS_ 2.979278E-07\n          62154._fuel_XS_ 1.935320E-06\n          63153._fuel_XS_ 5.206601E-06\n          63154._fuel_XS_ 6.662681E-07\n          63155._fuel_XS_ 2.851831E-07\n          63156._fuel_XS_ 4.669887E-07\n          64155._fuel_XS_ 5.962569E-10\n          64156._fuel_XS_ 1.009551E-06\n          64157._fuel_XS_ 9.145825E-09\n          64158._fuel_XS_ 5.483715E-07\nc  \nm13        6000._fuel_XS_ 1.182344e-02\n           8016._fuel_XS_ 3.547031e-02\n          90232._fuel_XS_ 3.377668E-12\n          90233._fuel_XS_ 2.448890E-17\n          91233._fuel_XS_ 5.920060E-13\n          92233._fuel_XS_ 1.462098E-11\n          92234._fuel_XS_ 2.454290E-08\n          92235._fuel_XS_ 2.167225E-03\n          92236._fuel_XS_ 4.681150E-04\n          92237._fuel_XS_ 3.689652E-06\n          92238._fuel_XS_ 1.792288E-02\n          92239._fuel_XS_ 1.495657E-07\n          93236._fuel_XS_ 2.250081E-12\n          93237._fuel_XS_ 2.646233E-05\n          93238._fuel_XS_ 2.774761E-07\n          93239._fuel_XS_ 2.130734E-05\n          94236._fuel_XS_ 5.636466E-12\n          94237._fuel_XS_ 1.758498E-12\n          94238._fuel_XS_ 4.367315E-06\n          94239._fuel_XS_ 2.908843E-04\n          94240._fuel_XS_ 8.973015E-05\n          94241._fuel_XS_ 7.834319E-05\n          94242._fuel_XS_ 1.862570E-05\n          94243._fuel_XS_ 1.927528E-08\n          94244._fuel_XS_ 1.179892E-09\n          95241._fuel_XS_ 3.882520E-07\n          95642._fuel_XS_ 2.472111E-08\n          95242._fuel_XS_ 3.709170E-09\n          95243._fuel_XS_ 2.291046E-06\n          95244._fuel_XS_ 3.763945E-10\n          96242._fuel_XS_ 1.243848E-07\n          96243._fuel_XS_ 1.560615E-09\n          96244._fuel_XS_ 3.989292E-07\n          96245._fuel_XS_ 1.690994E-08\n          96246._fuel_XS_ 7.294750E-10\n          96247._fuel_XS_ 4.768558E-12\n          96248._fuel_XS_ 1.691920E-13\n          97249._fuel_XS_ 1.334573E-15\n          98249._fuel_XS_ 3.437688E-17\n          98250._fuel_XS_ 4.957724E-16\n          35081._fuel_XS_ 5.092732E-06\n          36083._fuel_XS_ 1.087393E-05\n          36084._fuel_XS_ 2.548410E-05\n          36086._fuel_XS_ 4.553934E-05\n          37085._fuel_XS_ 2.296061E-05\n          37087._fuel_XS_ 5.808737E-05\n          38088._fuel_XS_ 8.355835E-05\n          38089._fuel_XS_ 4.030952E-05\n          38090._fuel_XS_ 1.281804E-04\n          39089._fuel_XS_ 6.719925E-05\n          39091._fuel_XS_ 5.677817E-05\n          40091._fuel_XS_ 7.767201E-05\n          40092._fuel_XS_ 1.798769E-04\n          40093._fuel_XS_ 1.513632E-04\n          40094._fuel_XS_ 1.507306E-04\n          40095._fuel_XS_ 7.095324E-05\n          40096._fuel_XS_ 1.546343E-04\n          42095._fuel_XS_ 5.067736E-05\n          42097._fuel_XS_ 1.501288E-04\n          42098._fuel_XS_ 1.495348E-04\n          42099._fuel_XS_ 3.417321E-06\n          42100._fuel_XS_ 1.630665E-04\n          43099._fuel_XS_ 1.376927E-04\n          44101._fuel_XS_ 1.302017E-04\n          44102._fuel_XS_ 1.211677E-04\n          44103._fuel_XS_ 3.307695E-05\n          44104._fuel_XS_ 6.868048E-05\n          44105._fuel_XS_ 9.637254E-08\n          44106._fuel_XS_ 2.647390E-05\n          45103._fuel_XS_ 5.129299E-05\n          45105._fuel_XS_ 5.901159E-07\n          46104._fuel_XS_ 1.303791E-05\n          46105._fuel_XS_ 3.504189E-05\n          46106._fuel_XS_ 1.527519E-05\n          46107._fuel_XS_ 1.983923E-05\n          46108._fuel_XS_ 1.321304E-05\n          46110._fuel_XS_ 3.990912E-06\n          47109._fuel_XS_ 6.917036E-06\n          48110._fuel_XS_ 1.802935E-06\n          48111._fuel_XS_ 1.899601E-06\n          48113._fuel_XS_ 1.744766E-08\n          48114._fuel_XS_ 1.504837E-06\n          49115._fuel_XS_ 2.811714E-07\n          52130._fuel_XS_ 1.000000E-30\n          53127._fuel_XS_ 4.716561E-06\n          53129._fuel_XS_ 2.110599E-05\n          54131._fuel_XS_ 5.765148E-05\n          54132._fuel_XS_ 1.265757E-04\n          54134._fuel_XS_ 1.966410E-04\n          54135._fuel_XS_ 4.840768E-08\n          54136._fuel_XS_ 3.248600E-04\n          55133._fuel_XS_ 1.472898E-04\n          55134._fuel_XS_ 1.386570E-05\n          55135._fuel_XS_ 1.759578E-05\n          55137._fuel_XS_ 1.571956E-04\n          56138._fuel_XS_ 1.702489E-04\n          56140._fuel_XS_ 1.617935E-05\n          57139._fuel_XS_ 1.601194E-04\n          58141._fuel_XS_ 3.857062E-05\n          58142._fuel_XS_ 1.488482E-04\n          58143._fuel_XS_ 1.531299E-06\n          59141._fuel_XS_ 1.054142E-04\n          59143._fuel_XS_ 1.548657E-05\n          60143._fuel_XS_ 1.053448E-04\n          60144._fuel_XS_ 4.541513E-05\n          60145._fuel_XS_ 8.836464E-05\n          60146._fuel_XS_ 8.207712E-05\n          60147._fuel_XS_ 4.852572E-06\n          60148._fuel_XS_ 4.633936E-05\n          60150._fuel_XS_ 1.826773E-05\n          61147._fuel_XS_ 3.014612E-05\n          61148._fuel_XS_ 8.658253E-07\n          61548._fuel_XS_ 2.123637E-07\n          61149._fuel_XS_ 7.482527E-07\n          62147._fuel_XS_ 1.698014E-06\n          62149._fuel_XS_ 4.153848E-07\n          62150._fuel_XS_ 3.132030E-05\n          62151._fuel_XS_ 1.881008E-06\n          62152._fuel_XS_ 1.236441E-05\n          62153._fuel_XS_ 4.217340E-07\n          62154._fuel_XS_ 3.557575E-06\n          63153._fuel_XS_ 1.042801E-05\n          63154._fuel_XS_ 1.707889E-06\n          63155._fuel_XS_ 5.794233E-07\n          63156._fuel_XS_ 9.442842E-07\n          64155._fuel_XS_ 1.212680E-09\n          64156._fuel_XS_ 2.929210E-06\n          64157._fuel_XS_ 1.433244E-08\n          64158._fuel_XS_ 1.259046E-06\nc  \nm14        6000._fuel_XS_ 1.182344e-02\n           8016._fuel_XS_ 3.547031e-02\n          90232._fuel_XS_ 5.892750E-12\n          90233._fuel_XS_ 4.267409E-17\n          91233._fuel_XS_ 1.137076E-12\n          92233._fuel_XS_ 1.523353E-11\n          92234._fuel_XS_ 3.563824E-08\n          92235._fuel_XS_ 1.551974E-03\n          92236._fuel_XS_ 5.594344E-04\n          92237._fuel_XS_ 4.543365E-06\n          92238._fuel_XS_ 1.752789E-02\n          92239._fuel_XS_ 1.484702E-07\n          93236._fuel_XS_ 4.188487E-12\n          93237._fuel_XS_ 4.463209E-05\n          93238._fuel_XS_ 4.753129E-07\n          93239._fuel_XS_ 2.120936E-05\n          94236._fuel_XS_ 1.280801E-11\n          94237._fuel_XS_ 4.262626E-12\n          94238._fuel_XS_ 1.057845E-05\n          94239._fuel_XS_ 2.990465E-04\n          94240._fuel_XS_ 1.035318E-04\n          94241._fuel_XS_ 1.087933E-04\n          94242._fuel_XS_ 3.973168E-05\n          94243._fuel_XS_ 4.081869E-08\n          94244._fuel_XS_ 4.134869E-09\n          95241._fuel_XS_ 7.039546E-07\n          95642._fuel_XS_ 5.552299E-08\n          95242._fuel_XS_ 6.766445E-09\n          95243._fuel_XS_ 6.970268E-06\n          95244._fuel_XS_ 1.163229E-09\n          96242._fuel_XS_ 3.392634E-07\n          96243._fuel_XS_ 5.917206E-09\n          96244._fuel_XS_ 1.798383E-06\n          96245._fuel_XS_ 1.003611E-07\n          96246._fuel_XS_ 6.321998E-09\n          96247._fuel_XS_ 5.631838E-11\n          96248._fuel_XS_ 2.849902E-12\n          97249._fuel_XS_ 2.708568E-14\n          98249._fuel_XS_ 8.959900E-16\n          98250._fuel_XS_ 1.221629E-14\n          35081._fuel_XS_ 6.537548E-06\n          36083._fuel_XS_ 1.306183E-05\n          36084._fuel_XS_ 3.347349E-05\n          36086._fuel_XS_ 5.792150E-05\n          37085._fuel_XS_ 2.929750E-05\n          37087._fuel_XS_ 7.379689E-05\n          38088._fuel_XS_ 1.060005E-04\n          38089._fuel_XS_ 3.622765E-05\n          38090._fuel_XS_ 1.622796E-04\n          39089._fuel_XS_ 1.000679E-04\n          39091._fuel_XS_ 5.293469E-05\n          40091._fuel_XS_ 1.183055E-04\n          40092._fuel_XS_ 2.254247E-04\n          40093._fuel_XS_ 1.946815E-04\n          40094._fuel_XS_ 1.954375E-04\n          40095._fuel_XS_ 6.963402E-05\n          40096._fuel_XS_ 2.014936E-04\n          42095._fuel_XS_ 9.202142E-05\n          42097._fuel_XS_ 1.968339E-04\n          42098._fuel_XS_ 1.971888E-04\n          42099._fuel_XS_ 2.927435E-06\n          42100._fuel_XS_ 2.155267E-04\n          43099._fuel_XS_ 1.766213E-04\n          44101._fuel_XS_ 1.709664E-04\n          44102._fuel_XS_ 1.652883E-04\n          44103._fuel_XS_ 3.333462E-05\n          44104._fuel_XS_ 9.808521E-05\n          44105._fuel_XS_ 9.773805E-08\n          44106._fuel_XS_ 3.962985E-05\n          45103._fuel_XS_ 7.268288E-05\n          45105._fuel_XS_ 5.731049E-07\n          46104._fuel_XS_ 2.822746E-05\n          46105._fuel_XS_ 5.158615E-05\n          46106._fuel_XS_ 2.585210E-05\n          46107._fuel_XS_ 3.192822E-05\n          46108._fuel_XS_ 2.175634E-05\n          46110._fuel_XS_ 6.585148E-06\n          47109._fuel_XS_ 1.081529E-05\n          48110._fuel_XS_ 4.111108E-06\n          48111._fuel_XS_ 3.168598E-06\n          48113._fuel_XS_ 1.824150E-08\n          48114._fuel_XS_ 2.237969E-06\n          49115._fuel_XS_ 3.265418E-07\n          52130._fuel_XS_ 1.000000E-30\n          53127._fuel_XS_ 6.656278E-06\n          53129._fuel_XS_ 2.870192E-05\n          54131._fuel_XS_ 7.144853E-05\n          54132._fuel_XS_ 1.775933E-04\n          54134._fuel_XS_ 2.592230E-04\n          54135._fuel_XS_ 4.129238E-08\n          54136._fuel_XS_ 4.306445E-04\n          55133._fuel_XS_ 1.891654E-04\n          55134._fuel_XS_ 2.542161E-05\n          55135._fuel_XS_ 2.406228E-05\n          55137._fuel_XS_ 2.075265E-04\n          56138._fuel_XS_ 2.230640E-04\n          56140._fuel_XS_ 1.384796E-05\n          57139._fuel_XS_ 2.094012E-04\n          58141._fuel_XS_ 3.446869E-05\n          58142._fuel_XS_ 1.951598E-04\n          58143._fuel_XS_ 1.271544E-06\n          59141._fuel_XS_ 1.530373E-04\n          59143._fuel_XS_ 1.307417E-05\n          60143._fuel_XS_ 1.323001E-04\n          60144._fuel_XS_ 8.265573E-05\n          60145._fuel_XS_ 1.113623E-04\n          60146._fuel_XS_ 1.113777E-04\n          60147._fuel_XS_ 4.165265E-06\n          60148._fuel_XS_ 6.126738E-05\n          60150._fuel_XS_ 2.471571E-05\n          61147._fuel_XS_ 3.319807E-05\n          61148._fuel_XS_ 9.745260E-07\n          61548._fuel_XS_ 2.353690E-07\n          61149._fuel_XS_ 7.184121E-07\n          62147._fuel_XS_ 2.756477E-06\n          62149._fuel_XS_ 3.974634E-07\n          62150._fuel_XS_ 4.325809E-05\n          62151._fuel_XS_ 2.045795E-06\n          62152._fuel_XS_ 1.502369E-05\n          62153._fuel_XS_ 4.908041E-07\n          62154._fuel_XS_ 5.315070E-06\n          63153._fuel_XS_ 1.572650E-05\n          63154._fuel_XS_ 2.941476E-06\n          63155._fuel_XS_ 9.401954E-07\n          63156._fuel_XS_ 1.544954E-06\n          64155._fuel_XS_ 2.055747E-09\n          64156._fuel_XS_ 6.374304E-06\n          64157._fuel_XS_ 2.084214E-08\n          64158._fuel_XS_ 2.300073E-06\nc\nm101 90232._fuel_XS_ 1.\nm102 90233._fuel_XS_ 1.\nm103 91233._fuel_XS_ 1.\nm104 92233._fuel_XS_ 1.\nm105 92234._fuel_XS_ 1.\nm106 92235._fuel_XS_ 1.\nm107 92236._fuel_XS_ 1.\nm108 92237._fuel_XS_ 1.\nm109 92238._fuel_XS_ 1.\nm110 92239._fuel_XS_ 1.\nm111 93236._fuel_XS_ 1.\nm112 93237._fuel_XS_ 1.\nm113 93238._fuel_XS_ 1.\nm114 93239._fuel_XS_ 1.\nm115 94236._fuel_XS_ 1.\nm116 94237._fuel_XS_ 1.\nm117 94238._fuel_XS_ 1.\nm118 94239._fuel_XS_ 1.\nm119 94240._fuel_XS_ 1.\nm120 94241._fuel_XS_ 1.\nm121 94242._fuel_XS_ 1.\nm122 94243._fuel_XS_ 1.\nm123 94244._fuel_XS_ 1.\nm124 95241._fuel_XS_ 1.\nm125 95642._fuel_XS_ 1.\nm126 95242._fuel_XS_ 1.\nm127 95243._fuel_XS_ 1.\nm128 95244._fuel_XS_ 1.\nm129 96242._fuel_XS_ 1.\nm130 96243._fuel_XS_ 1.\nm131 96244._fuel_XS_ 1.\nm132 96245._fuel_XS_ 1.\nm133 96246._fuel_XS_ 1.\nm134 96247._fuel_XS_ 1.\nm135 96248._fuel_XS_ 1.\nm136 97249._fuel_XS_ 1.\nm137 98249._fuel_XS_ 1.\nm138 98250._fuel_XS_ 1.\nm139 35081._fuel_XS_ 1.\nm140 36083._fuel_XS_ 1.\nm141 36084._fuel_XS_ 1.\nm142 36086._fuel_XS_ 1.\nm143 37085._fuel_XS_ 1.\nm144 37087._fuel_XS_ 1.\nm145 38088._fuel_XS_ 1.\nm146 38089._fuel_XS_ 1.\nm147 38090._fuel_XS_ 1.\nm148 39089._fuel_XS_ 1.\nm149 39091._fuel_XS_ 1.\nm150 40091._fuel_XS_ 1.\nm151 40092._fuel_XS_ 1.\nm152 40093._fuel_XS_ 1.\nm153 40094._fuel_XS_ 1.\nm154 40095._fuel_XS_ 1.\nm155 40096._fuel_XS_ 1.\nm156 42095._fuel_XS_ 1.\nm157 42097._fuel_XS_ 1.\nm158 42098._fuel_XS_ 1.\nm159 42099._fuel_XS_ 1.\nm160 42100._fuel_XS_ 1.\nm161 43099._fuel_XS_ 1.\nm162 44101._fuel_XS_ 1.\nm163 44102._fuel_XS_ 1.\nm164 44103._fuel_XS_ 1.\nm165 44104._fuel_XS_ 1.\nm166 44105._fuel_XS_ 1.\nm167 44106._fuel_XS_ 1.\nm168 45103._fuel_XS_ 1.\nm169 45105._fuel_XS_ 1.\nm170 46104._fuel_XS_ 1.\nm171 46105._fuel_XS_ 1.\nm172 46106._fuel_XS_ 1.\nm173 46107._fuel_XS_ 1.\nm174 46108._fuel_XS_ 1.\nm175 46110._fuel_XS_ 1.\nm176 47109._fuel_XS_ 1.\nm177 48110._fuel_XS_ 1.\nm178 48111._fuel_XS_ 1.\nm179 48113._fuel_XS_ 1.\nm180 48114._fuel_XS_ 1.\nm181 49115._fuel_XS_ 1.\nm182 52130._fuel_XS_ 1.\nm183 53127._fuel_XS_ 1.\nm184 53129._fuel_XS_ 1.\nm185 54131._fuel_XS_ 1.\nm186 54132._fuel_XS_ 1.\nm187 54134._fuel_XS_ 1.\nm188 54135._fuel_XS_ 1.\nm189 54136._fuel_XS_ 1.\nm190 55133._fuel_XS_ 1.\nm191 55134._fuel_XS_ 1.\nm192 55135._fuel_XS_ 1.\nm193 55137._fuel_XS_ 1.\nm194 56138._fuel_XS_ 1.\nm195 56140._fuel_XS_ 1.\nm196 57139._fuel_XS_ 1.\nm197 58141._fuel_XS_ 1.\nm198 58142._fuel_XS_ 1.\nm199 58143._fuel_XS_ 1.\nm200 59141._fuel_XS_ 1.\nm201 59143._fuel_XS_ 1.\nm202 60143._fuel_XS_ 1.\nm203 60144._fuel_XS_ 1.\nm204 60145._fuel_XS_ 1.\nm205 60146._fuel_XS_ 1.\nm206 60147._fuel_XS_ 1.\nm207 60148._fuel_XS_ 1.\nm208 60150._fuel_XS_ 1.\nm209 61147._fuel_XS_ 1.\nm210 61148._fuel_XS_ 1.\nm211 61548._fuel_XS_ 1.\nm212 61149._fuel_XS_ 1.\nm213 62147._fuel_XS_ 1.\nm214 62149._fuel_XS_ 1.\nm215 62150._fuel_XS_ 1.\nm216 62151._fuel_XS_ 1.\nm217 62152._fuel_XS_ 1.\nm218 62153._fuel_XS_ 1.\nm219 62154._fuel_XS_ 1.\nm220 63153._fuel_XS_ 1.\nm221 63154._fuel_XS_ 1.\nm222 63155._fuel_XS_ 1.\nm223 63156._fuel_XS_ 1.\nm224 64155._fuel_XS_ 1.\nm225 64156._fuel_XS_ 1.\nm226 64157._fuel_XS_ 1.\nm227 64158._fuel_XS_ 1.\nc\nm1  6000._core_XS_  1.\nmt1 grph.64t\nm2  6000._shell_XS_ 1.\nmt2 grph.64t \nc   Coatings and matrix mix\nm19    6000._fuel_XS_  _seed_matrix_carbon_\n      14028._fuel_XS_  _seed_matrix_silicon_\nmt19    grph.64t\nc   Salt in the Core - 2LiF-BeF2 - (99.995w% Li-6)\nm30    3006._coolant_XS_  1.166657e-04\n       3007._coolant_XS_  1.999883\n       4009._coolant_XS_  1.\n       9019._coolant_XS_  4.\nc Steel for Reactor Pressure Vessel and Core Barrel\nm31       28058.71c 0.0438374141756552\n          28060.71c 0.0163231374502\n          28061.71c 0.00069798731409836\n          28062.71c 0.00218910073380645\n          28064.71c 0.000540385299125\n          24050.71c 0.00032455797836\n          24052.71c 0.00601805533522308\n          24053.71c 0.000669523710271698\n          24054.71c 0.000163572351962963\n          42092.71c 0.00154914850330435\n          42094.71c 0.000945063210638298\n          42095.71c 0.00160940903949474\n          42096.71c 0.001668675123\n          42097.71c 0.000945537190515464\n          42098.71c 0.00236471183191837\n          42100.71c 0.00092485418328\n          26054.71c 0.000286535681481481\n          26056.71c 0.00436938684142857\n          26057.71c 0.000102965635087719\n          26058.71c 1.28787737931034E-05\nc    begin_mocup_reaction_rate_tallies\nc    time dependent reaction rates\nfc104  Reaction rates\nf104:n   \n       10 15 20 25   \nfm104  (1)\n       (1 101 (16) (17) (-6) (102) )\n       (1 102 (16) (17) (-6) (102) )\n       (1 103 (16) (17) (-6) (102) )\n       (1 104 (16) (17) (-6) (102) )\n       (1 105 (16) (17) (-6) (102) )\n       (1 106 (16) (17) (-6) (102) )\n       (1 107 (16) (17) (-6) (102) )\n       (1 108 (16) (17) (-6) (102) )\n       (1 109 (16) (17) (-6) (102) )\n       (1 110 (16) (17) (-6) (102) )\n       (1 111 (16) (17) (-6) (102) )\n       (1 112 (16) (17) (-6) (102) )\n       (1 113 (16) (17) (-6) (102) )\n       (1 114 (16) (17) (-6) (102) )\n       (1 115 (16) (17) (-6) (102) )\n       (1 116 (16) (17) (-6) (102) )\n       (1 117 (16) (17) (-6) (102) )\n       (1 118 (16) (17) (-6) (102) )\n       (1 119 (16) (17) (-6) (102) )\n       (1 120 (16) (17) (-6) (102) )\n       (1 121 (16) (17) (-6) (102) )\n       (1 122 (16) (17) (-6) (102) )\n       (1 123 (16) (17) (-6) (102) )\n       (1 124 (16) (17) (-6) (102) )\n       (1 125 (16) (17) (-6) (102) )\n       (1 126 (16) (17) (-6) (102) )\n       (1 127 (16) (17) (-6) (102) )\n       (1 128 (16) (17) (-6) (102) )\n       (1 129 (16) (17) (-6) (102) )\n       (1 130 (16) (17) (-6) (102) )\n       (1 131 (16) (17) (-6) (102) )\n       (1 132 (16) (17) (-6) (102) )\n       (1 133 (16) (17) (-6) (102) )\n       (1 134 (16) (17) (-6) (102) )\n       (1 135 (16) (17) (-6) (102) )\n       (1 136 (16) (17) (-6) (102) )\n       (1 137 (16) (17) (-6) (102) )\n       (1 138 (16) (17) (-6) (102) )\n       (1 139 (16) (17)      (102) )\n       (1 140 (16) (17)      (102) )\n       (1 141 (16) (17)      (102) )\n       (1 142 (16) (17)      (102) )\n       (1 143 (16) (17)      (102) )\n       (1 144 (16) (17)      (102) )\n       (1 145 (16) (17)      (102) )\n       (1 146 (16) (17)      (102) )\n       (1 147 (16) (17)      (102) )\n       (1 148 (16) (17)      (102) )\n       (1 149 (16) (17)      (102) )\n       (1 150 (16) (17)      (102) )\n       (1 151 (16) (17)      (102) )\n       (1 152 (16) (17)      (102) )\n       (1 153 (16) (17)      (102) )\n       (1 154 (16) (17)      (102) )\n       (1 155 (16) (17)      (102) )\n       (1 156 (16) (17)      (102) )\n       (1 157 (16) (17)      (102) )\n       (1 158 (16) (17)      (102) )\n       (1 159 (16) (17)      (102) )\n       (1 160 (16) (17)      (102) )\n       (1 161 (16) (17)      (102) )\n       (1 162 (16) (17)      (102) )\n       (1 163 (16) (17)      (102) )\n       (1 164 (16) (17)      (102) )\n       (1 165 (16) (17)      (102) )\n       (1 166 (16) (17)      (102) )\n       (1 167 (16) (17)      (102) )\n       (1 168 (16) (17)      (102) )\n       (1 169 (16) (17)      (102) )\n       (1 170 (16) (17)      (102) )\n       (1 171 (16) (17)      (102) )\n       (1 172 (16) (17)      (102) )\n       (1 173 (16) (17)      (102) )\n       (1 174 (16) (17)      (102) )\n       (1 175 (16) (17)      (102) )\n       (1 176 (16) (17)      (102) )\n       (1 177 (16) (17)      (102) )\n       (1 178 (16) (17)      (102) )\n       (1 179 (16) (17)      (102) )\n       (1 180 (16) (17)      (102) )\n       (1 181 (16) (17)      (102) )\n       (1 182 (16) (17)      (102) )\n       (1 183 (16) (17)      (102) )\n       (1 184 (16) (17)      (102) )\n       (1 185 (16) (17)      (102) )\n       (1 186 (16) (17)      (102) )\n       (1 187 (16) (17)      (102) )\n       (1 188 (16) (17)      (102) )\n       (1 189 (16) (17)      (102) )\n       (1 190 (16) (17)      (102) )\n       (1 191 (16) (17)      (102) )\n       (1 192 (16) (17)      (102) )\n       (1 193 (16) (17)      (102) )\n       (1 194 (16) (17)      (102) )\n       (1 195 (16) (17)      (102) )\n       (1 196 (16) (17)      (102) )\n       (1 197 (16) (17)      (102) )\n       (1 198 (16) (17)      (102) )\n       (1 199 (16) (17)      (102) )\n       (1 200 (16) (17)      (102) )\n       (1 201 (16) (17)      (102) )\n       (1 202 (16) (17)      (102) )\n       (1 203 (16) (17)      (102) )\n       (1 204 (16) (17)      (102) )\n       (1 205 (16) (17)      (102) )\n       (1 206 (16) (17)      (102) )\n       (1 207 (16) (17)      (102) )\n       (1 208 (16) (17)      (102) )\n       (1 209 (16) (17)      (102) )\n       (1 210 (16) (17)      (102) )\n       (1 211 (16) (17)      (102) )\n       (1 212 (16) (17)      (102) )\n       (1 213 (16) (17)      (102) )\n       (1 214 (16) (17)      (102) )\n       (1 215 (16) (17)      (102) )\n       (1 216 (16) (17)      (102) )\n       (1 217 (16) (17)      (102) )\n       (1 218 (16) (17)      (102) )\n       (1 219 (16) (17)      (102) )\n       (1 220 (16) (17)      (102) )\n       (1 221 (16) (17)      (102) )\n       (1 222 (16) (17)      (102) )\n       (1 223 (16) (17)      (102) )\n       (1 224 (16) (17)      (102) )\n       (1 225 (16) (17)      (102) )\n       (1 226 (16) (17)      (102) )\n       (1 227 (16) (17)      (102) )\nfc2 tally for baseline flux out of reactor\n f2:n      803\nfs2      -923\nsd2       _area2_ 1\nfc22 tally of fast fluence to outer graphite reflector\n f22:n 923\nfs22 803 -802\n e22 .1 20\nsd22 1 1 1 \nfc6 tally for power fraction in converging section\n f6:n     ((10 15 20 25) < 931) \n           (10 15 20 25)\nsd6 1 1 \nkcode  10000  1.0 10 110\nprdmp  10000 10000 10000\n#print\nmode n\nksrc 0 0 50\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n'

        self.wall_effects_radius = self.characteristic_radius - \
            self.Pebble.r_shell*self.wall_effects_normalized
        self.shutdown_rod_r = self.characteristic_radius * \
            math.tan(math.pi/(2*self.shutdown_rods))/(1 - 2*math.tan(math.pi/(2*self.shutdown_rods)))

        self.dr_graphite_reflector = 4*self.shutdown_rod_r
        self.graphite_reflector_radius = self.characteristic_radius + \
            self.dr_graphite_reflector
        self.barrel_radius = self.graphite_reflector_radius + self.dr_barrel

        # the flow area of coolant through the core
        flow_area = math.pi*(math.pow(self.characteristic_radius,
                                      2.) - math.pow(self.wall_effects_radius,
                                                     2.))*self.high_porosity
        flow_area += math.pi * \
            math.pow(self.wall_effects_radius, 2.)*self.low_porosity

        self.down_comer_radius = math.pow(
            ((flow_area + math.pi * math.pow(self.barrel_radius, 2.)) / math.pi),
            0.5)
        self.pressure_vessel_radius = self.down_comer_radius + \
            self.dr_pressure_vessel
        self.dz_buffer = active_volume/(1. - self.defueling_chute_volume_normalized) * \
            self.buffer_volume_normalized/(math.pi*math.pow(self.characteristic_radius, 2.))

        self.defueling_chute_radius = self.defueling_radius_normalized * \
            self.Pebble.r_shell
        self.dz_converging = (
            self.characteristic_radius - self.defueling_chute_radius)*math.tan(
            (self.converging_angle*math.pi/180.))
        self.dz_defueling_chute = self.defueling_chute_volume_normalized*active_volume / \
            (1. - self.defueling_chute_volume_normalized)/(math.pi*math.pow(self.defueling_chute_radius, 2.))

        dz_plenum = self.down_comer_radius*2/self.lower_plenum_aspect_ratio

        conical_intersect = self.dz_buffer + self.active_height + \
            self.characteristic_radius*math.tan(self.converging_angle*math.pi/180.)
        tan_squared = math.pow(
            math.tan(2*math.pi-(math.pi/180.)*self.converging_angle), 2.)

        volume = (active_volume / (1. - self.defueling_chute_volume_normalized) -
                  math.pi * math.pow(self.wall_effects_radius, 2.) *
                  self.active_height) * (1. - self.high_porosity)
        volume += math.pi * \
            math.pow(self.wall_effects_radius, 2.)*self.active_height*(1.-self.low_porosity)

        self.volume = volume*(self.Pebble.dv_active/self.Pebble.V) * \
            (self.Pebble.TRISO.dv_kernel/self.Pebble.TRISO.V)/4.

        CHM, rho = self.Pebble.CHM()
        self.Pebble.TRISO.matrix()
        # generate outer reflector cells

        theta = 2.*math.pi/self.shutdown_rods
        radius = self.characteristic_radius + 2.*self.shutdown_rod_r
        _SR_Channel_Surfaces_ = 'c Surfaces for Shutdown Channels'
        _Outer_Reflector_ = '800 1 -1.74 923 -904 -805 812 '
        _Shutdown_Channels_ = 'c Shutdown Channels'

        for i in range(self.shutdown_rods):
            X = radius*math.sin(i*theta)
            Y = radius*math.cos(i*theta)
            surface_number = 8000 + i
            _SR_Channel_Surfaces_ += '\n%d c/z %1.5E %1.5E %1.5E' % (
                surface_number, X, Y, self.shutdown_rod_r)
            _Outer_Reflector_ += '\n          %d' % (8000 + i)
            _Shutdown_Channels_ += '\n%d 30 -_coolant_density_ -805 812 -%d imp:n=1 tmp=%1.5E' % (
                (8000 + i), (8000 + i), (self.T_bulk*8.6173324e-11))
        _Outer_Reflector_ += '\n          imp:n=1'

        # _plenum1_
        _plenum1_ = '%1.5E %1.5E %1.5E' % (
            (1. / math.pow(self.pressure_vessel_radius, 2.)),
            (1. / math.pow(self.pressure_vessel_radius, 2.)),
            (1. / math.pow((dz_plenum + self.dr_pressure_vessel), 2.)))
        # _plenum2_
        _plenum2_ = '%1.5E %1.5E %1.5E' % (
            (1. / math.pow(self.down_comer_radius, 2.)),
            (1. / math.pow(self.down_comer_radius, 2.)),
            (1. / math.pow(dz_plenum, 2.)))

        _area2_ = '%1.5E' % (
            math.pi*self.characteristic_radius*self.active_height)

        if self.Pebble.T_core < (750):
            _core_XS_ = '71c'
        elif self.Pebble.T_core < (1050) and self.Pebble.T_core > (750):
            _core_XS_ = '72c'
        else:
            _core_XS_ = '73c'

        if self.Pebble.T_kernel < (750):
            _fuel_XS_ = '71c'
        elif self.Pebble.T_kernel < (1050) and self.Pebble.T_kernel > (750):
            _fuel_XS_ = '72c'
        else:
            _fuel_XS_ = '73c'

        if self.Pebble.T_shell < (750):
            _shell_XS_ = '71c'
        elif self.Pebble.T_shell < (1050) and self.Pebble.T_shell > (750):
            _shell_XS_ = '72c'
        else:
            _shell_XS_ = '73c'

        if self.T_bulk < (750):
            _coolant_XS_ = '71c'
        elif self.T_bulk < (1050) and self.T_bulk > (750):
            _coolant_XS_ = '72c'
        else:
            _coolant_XS_ = '73c'

        _seed_title_ = 'c Test Reactor\nc TRISO Fuel Design - kernel diameter %4.1f (um) | C/HM %4.1f ' % (
            self.Pebble.TRISO.r_kernel*1e4, CHM)
        _fuel_temperature_ = '%1.5E' % (self.Pebble.T_kernel*8.6173324e-11)
        _seed_matrix_rho_ = '%1.5E' % self.Pebble.TRISO.mat_hmatrix.mass()
        _seed_kernel_volume_ = '%1.5E' % self.volume
        _coolant_temperature_ = '%1.5E' % (self.T_bulk*8.6173324e-11)
        _seed_core_rho_ = '%1.5E' % self.Pebble.mat_core.mass()
        _seed_core_temperature_ = '%1.5E' % (self.Pebble.T_core*8.6173324e-11)
        _seed_shell_rho_ = '%1.5E' % self.Pebble.mat_shell.mass()
        _seed_shell_temperature_ = '%1.5E' % (self.Pebble.T_shell*8.6173324e-11)
        _coolant_density_ = '%1.5E' % (2.28-.000488*(self.T_bulk-273.15))
        _seed_kernel_r_ = '%1.5E' % self.Pebble.TRISO.r_kernel
        _seed_TRISO_hp_ = '%1.5E' % self.Pebble.TRISO.hpitch
        _seed_pebble_hp1_ = '%1.5E' % (
            math.pow(4*self.Pebble.V/(1.-self.low_porosity), (1./3.))/2.)
        _seed_pebble_hp2_ = '%1.5E' % (
            math.pow(4*self.Pebble.V/(1.-self.high_porosity), (1./3.))/2.)
        _seed_core_r_ = '%1.5E' % self.Pebble.r_core
        _seed_active_r_ = '%1.5E' % self.Pebble.r_active
        _seed_shell_r_ = '%1.5E' % self.Pebble.r_shell
        _OGRR_ = '%1.5E' % self.graphite_reflector_radius
        _core_barrel_r_ = '%1.5E' % self.barrel_radius
        _down_comer_r_ = '%1.5E' % self.down_comer_radius
        _pressure_vessel_r_ = '%1.5E' % self.pressure_vessel_radius
        _boundary_layer_R_ = '%1.5E' % self.wall_effects_radius
        _active_R_ = '%1.5E' % self.characteristic_radius
        _contraction_b_ = '%1.5E' % conical_intersect
        _contraction_tan_ = '%1.5E' % tan_squared
        _exit_R_ = '%1.5E' % self.defueling_chute_radius
        _free_surface_ = '%1.5E' % (self.dz_buffer)
        _active_Z_ = '%1.5E' % (self.dz_buffer + self.active_height)
        _contraction_Z_ = '%1.5E' % (
            self.dz_buffer + self.active_height + self.dz_converging)
        _exit_Z_ = '%1.5E' % (
            self.dz_buffer + self.active_height + self.dz_converging +
            self.dz_defueling_chute)
        _seed_matrix_carbon_ = '%1.5E' % (
            self.Pebble.TRISO.mat_hmatrix.comp['60120']/12.)
        _seed_matrix_silicon_ = '%1.5E' % (
            self.Pebble.TRISO.mat_hmatrix.comp['140280']/28.)

        skeleton = skeleton.replace('_Outer_Reflector_', _Outer_Reflector_)
        skeleton = skeleton.replace('_Shutdown_Channels_', _Shutdown_Channels_)
        skeleton = skeleton.replace(
            '_SR_Channel_Surfaces_',
            _SR_Channel_Surfaces_)
        skeleton = skeleton.replace('_seed_title_', _seed_title_)
        skeleton = skeleton.replace('_fuel_temperature_', _fuel_temperature_)
        skeleton = skeleton.replace('_seed_matrix_rho_', _seed_matrix_rho_)
        skeleton = skeleton.replace(
            '_coolant_temperature_',
            _coolant_temperature_)
        skeleton = skeleton.replace(
            '_seed_kernel_volume_',
            _seed_kernel_volume_)
        skeleton = skeleton.replace('_seed_core_rho_', _seed_core_rho_)
        skeleton = skeleton.replace(
            '_seed_core_temperature_',
            _seed_core_temperature_)
        skeleton = skeleton.replace('_seed_shell_rho_', _seed_shell_rho_)
        skeleton = skeleton.replace(
            '_seed_shell_temperature_',
            _seed_shell_temperature_)

        skeleton = skeleton.replace('_coolant_density_', _coolant_density_)
        skeleton = skeleton.replace('_seed_kernel_r_', _seed_kernel_r_)
        skeleton = skeleton.replace('_seed_TRISO_hp_', _seed_TRISO_hp_)
        skeleton = skeleton.replace('_seed_pebble_hp1_', _seed_pebble_hp1_)
        skeleton = skeleton.replace('_seed_pebble_hp2_', _seed_pebble_hp2_)
        skeleton = skeleton.replace('_seed_core_r_', _seed_core_r_)
        skeleton = skeleton.replace('_seed_active_r_', _seed_active_r_)
        skeleton = skeleton.replace('_seed_shell_r_', _seed_shell_r_)
        skeleton = skeleton.replace('_OGRR_', _OGRR_)
        skeleton = skeleton.replace('_core_barrel_r_', _core_barrel_r_)

        skeleton = skeleton.replace('_down_comer_r_', _down_comer_r_)
        skeleton = skeleton.replace('_pressure_vessel_r_', _pressure_vessel_r_)
        skeleton = skeleton.replace('_boundary_layer_R_', _boundary_layer_R_)
        skeleton = skeleton.replace('_active_R_', _active_R_)
        skeleton = skeleton.replace('_contraction_b_', _contraction_b_)
        skeleton = skeleton.replace('_contraction_tan_', _contraction_tan_)
        skeleton = skeleton.replace('_exit_R_', _exit_R_)
        skeleton = skeleton.replace('_plenum1_', _plenum1_)
        skeleton = skeleton.replace('_plenum2_', _plenum2_)
        skeleton = skeleton.replace('_free_surface_', _free_surface_)

        skeleton = skeleton.replace('_active_Z_', _active_Z_)
        skeleton = skeleton.replace('_contraction_Z_', _contraction_Z_)
        skeleton = skeleton.replace('_exit_Z_', _exit_Z_)
        skeleton = skeleton.replace(
            '_seed_matrix_carbon_',
            _seed_matrix_carbon_)
        skeleton = skeleton.replace(
            '_seed_matrix_silicon_',
            _seed_matrix_silicon_)
        skeleton = skeleton.replace('_core_XS_', _core_XS_)
        skeleton = skeleton.replace('_fuel_XS_', _fuel_XS_)
        skeleton = skeleton.replace('_shell_XS_', _shell_XS_)
        skeleton = skeleton.replace('_coolant_XS_', _coolant_XS_)
        skeleton = skeleton.replace('_area2_', _area2_)

        self.input = skeleton
        open('inp.1', 'w').write(skeleton)

    def BEAU(self):

        # This module sets up a Equilibrium Depletion analysis directory

        import os

        self.inputMCNP5()
        dir = '../Equilibrium_Test_%s' % self.title

        # generate a BEAU executable for the specific test reactor design
        BEAU = "#!/usr/bin/env python\n\nimport BEAU\nimport mocup\n\n# -----------------------------------------------------------------------------\n# DEFINE EQUILIBRIUM DEPLETION ANALYSIS PROBLEM\n# -----------------------------------------------------------------------------\n\n# initialize Equilibrium Depletion\nEq_Cycle = BEAU.equilibrium_cycle()\n\n# DEFINE POWER OF SYSTEM IN MWth\nEq_Cycle.power = _power_\n\t# 8 1/2 pebbles with a volume of 1.17810E+01 cc with a power desnity of 30 MW/m3\n\n# DEFINE LOADING PATTER\n\t# 'Y':'X' EOEC material from X is advanced to cell Y\nEq_Cycle.loading_pattern = {\n\t\t\t\t\t\t\t'15':'10',\n\t\t\t\t\t\t\t'20':'15',\n\t\t\t\t\t\t\t'25':'20',\n\t\t\t\t\t\t\t}\n\n# DEFINE THE BURNUPS OF ALL FUEL PROGRESSIONS\nEq_Cycle.burnup = {'seed': 1.5E+5}\n\n# DEFINE THE MAKE UP MATERIAL FOR EACH PROGRESSION\n\n# initiate material object\nmat = mocup.material()\n# define its composition vector\n_material_\n\n# link material object to the progression key\nEq_Cycle.makeup = {'seed':mat}\n\n# INITIATE A PROGRESSION\nprogression = Eq_Cycle.fuel_progression()\n\n# DEFINE THE ORDER IN WHICH FUEL IS ADVANCED\nprogression.cells = ['10','15','20','25']\n\n# DEFINE THE MAKEUP FUEL VECTOR FOR THE PROGRESSION\nprogression.feed = Eq_Cycle.makeup['seed']\n\n# DEFINE THE STANDARD ORIGEN LIBRARY TO BE USED\nprogression.library = 'pwru50'\n\n# LINK PROGRESSION OBJECTS TO KEYS\nEq_Cycle.progressions = {'seed':progression,}\n\n# -----------------------------------------------------------------------------\n# INITIATE EQUILIBRIUM DEPLETION ANALYSIS\n# -----------------------------------------------------------------------------\nEq_Cycle.dBU1 = 1000\nEq_Cycle.dBU2 = 10000\nEq_Cycle.nBU1 = 3\n\nEq_Cycle.search_keff(1.0E+05,1.5e+05, 1.000)\n"
        _power_ = '%1.5E' % self.power
        _material_ = ''
        for nucl in list(self.Pebble.TRISO.mat_kernel.comp.keys()):
            if float(nucl[:-4]) > 80.:
                # this nuclide is a heavy metal nuclide
                _material_ += "\nmat.comp['%s'] = %1.5E" % (nucl,
                                                            self.volume*self.Pebble.TRISO.mat_kernel.comp[nucl])
        BEAU = BEAU.replace('_power_', _power_)
        BEAU = BEAU.replace('_material_', _material_)

        # copy Template directory to main directory

        copy = '\\rm -r %s' % dir
        # print(copy)
        os.system(copy)

        copy = '\cp -r Equilibrium_BEAU %s' % dir
        # print(copy)
        os.system(copy)

        # move input deck and executable to the equilibrium depletion analysis
        # directory

        inp_loc = dir + '/inp.1'
        BEAU_loc = dir + '/BEAU7.py'

        copy = '\cp inp.1 %s' % inp_loc
        # print(copy)
        os.system(copy)

        open(BEAU_loc, 'w').write(BEAU)

        chmod = 'chmod +x %s' % BEAU_loc
        # print(chmod)
        os.system(chmod)

        # comment out AUSCE

        BEAU_loc = dir + '/BEAU.py'
        BEAU = open(BEAU_loc).read()
        BEAU = BEAU.replace('self.AUSCE(', '#self.AUSCE(')
        open(BEAU_loc, 'w').write(BEAU)

        mcnp = "#! /bin/sh\n#\n#$ -N Equilibrium_Depletion\n#$ -cwd\n#$ -pe ompi 8\n#$ -S /bin/bash\n#$ -q x.q\n#$ -V\n\nmpiexec -x LD_LIBRARY_PATH  mcnp5.mpi i=inp.1 o=outp.1 mc=mctal.1 runtpe=runtp1\nmv srctp source\nrm outp.1 mctal.1 tmp runtp1\nsed s/'ksrc 0 0 50'/'c'/ inp.1 >> tmp\nmv tmp inp.1\n\n./BEAU7.py\n\n"

        mcnp_loc = dir + '/mcnp.sh'
        open(mcnp_loc, 'w').write(mcnp)

        chmod = 'chmod +x %s' % (mcnp_loc)
        # print(chmod)
        os.system(chmod)

