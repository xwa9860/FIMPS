#!/usr/bin/python


class ANNULAR_Graphite:

    def __init__(self):
        # initialized the seed fuel design
        self.seed = Pebble()

        # set the seed TRISO fuel design to baseline
        self.seed.TRISO.dr_kernel = .025
        self.seed.TRISO.setPF(.11121)

        # set the seed Pebble fuel design to baseline
        self.seed.dr_core = .577
        self.seed.dr_shell = .1
        self.seed.volume()
        self.seed.mat_core.comp['60120'] = (.5/12.)
        self.seed.TRISO.dr_kernel = .02
        self.seed.generate(300., 1.745)

        # initialize the blanket fuel design
        self.blanket = Pebble()
        self.blanket.mat_shell.comp['60120'] = 1.74/12.
        self.blanket.volume()

        # initialize shield
        self.absorber = shield()

        self.min_annular_layer_normalized = 4.
        self.inner_reflector_radius = 90.
        self.characteristic_radius = 240.
        self.dz_active = 300.
        self.dr_outer_reflector = 50.
        self.inlet_annular_layer_normalized = 23.3333
        self.outlet_annular_layer_normalized = 10.
        self.expansion_angle = 60.
        self.converging_angle = 45.
        self.z_free_surface = 0.
        self.dz_entrance = 150.
        self.dz_defueling = 150.
        self.plenum_aspect_ratio = 1./4.
        self.superficial_angle = 60.
        self.CTE = 1.e-04

        self.inner_shutdown_rods = 12
        # radius of inner shutdown rods is based on maximum radius with robust constructruction (one radius between perimeter of channel and any other boundaries, ie other channels or the active annulus)
        #self.outer_shutdown_rods = 32
        self.outer_shutdown_rods = 0
        # radius of inner shutdown rods is based on maximum radius with robust constructruction (one radius between perimeter of channel and any other boundaries, ie other channels or the active annulus)
        #self.shutdown_rod_efficiency = 0.779914334
        self.shutdown_rod_efficiency = 1e-2

        self.dr_core_barrel = .3
        self.dr_pressure_vessel = .3

        self.power_density = 12.78689745*1.266072082
        self.outlet_temperature = (700 + 273.15)
        self.temperature_rise = 100.

        self.dr_blanket_entrance = self.min_annular_layer_normalized * \
            self.blanket.r_shell*2.
        self.assumed_blanket_power = 0

    def no_shield(self):
        # set the pitch of the absorbers so larger that they will likely not be present in the calculation
        # there will be only one absorber pin at the center of the system, which
        # is not where the shield looks into the infinite array of absorber rods
        import math
        PF = math.pow(self.absorber.absorber_radius, 2.) / \
            math.pow(self.characteristic_radius, 2.)
        self.absorber.setPF(PF)

    def min_blanket(self):
        self.dr_blanket_entrance = self.min_annular_layer_normalized * \
            self.blanket.r_shell*2.

    def generate(self):
        # This method generates the geometry of the reactor

        import math

        # calculate the seed-blanket interfaces in each of the cylindrical
        # regions

        self.blanket_entrance_radius = self.inner_reflector_radius + \
            self.inlet_annular_layer_normalized*self.blanket.r_shell*2.
        self.seed_entrance_radius = self.blanket_entrance_radius - \
            self.dr_blanket_entrance

        seed_fraction = (math.pow(self.seed_entrance_radius,
                                  2.) - math.pow(self.inner_reflector_radius,
                                                 2.))/(math.pow(self.blanket_entrance_radius,
                                                                2.) - math.pow(self.inner_reflector_radius,
                                                                               2.))
        blanket_fraction = 1 - seed_fraction

        reflector_area = math.pi*math.pow(self.inner_reflector_radius, 2.)
        active_area = math.pi * \
            math.pow(self.characteristic_radius, 2.) - reflector_area

        self.seed_active_radius = math.pow(
            (reflector_area + seed_fraction*active_area)/math.pi,
            0.5)
        self.blanket_active_radius = self.characteristic_radius

        self.blanket_defueling_radius = self.inner_reflector_radius + \
            self.seed.r_shell*2*self.outlet_annular_layer_normalized
        defueling_area = math.pi * \
            math.pow(self.blanket_defueling_radius, 2.) - reflector_area
        self.seed_defueling_radius = math.pow(
            (reflector_area +
             seed_fraction *
             defueling_area) /
            math.pi,
            0.5)

        # calculate the angles and z-intercepts of the cones

        self.z_entrance = self.z_free_surface + self.dz_entrance

        # expansion region

        self.dz_expansion = math.tan(
            self.expansion_angle*math.pi/180.)*(self.blanket_active_radius - self.blanket_entrance_radius)
        m = (self.dz_expansion) / \
            (self.seed_active_radius - self.seed_entrance_radius)
        self.seed_expansion_intercept = self.z_entrance - \
            m*self.seed_entrance_radius
        adjacent = self.z_entrance - self.seed_expansion_intercept
        seed_expansion_angle = 2.*math.atan(self.seed_entrance_radius/adjacent)

        m = math.tan(self.expansion_angle*math.pi/180.)
        self.blanket_expansion_intercept = self.z_entrance - \
            m*self.blanket_entrance_radius
        adjacent = self.z_entrance - self.blanket_expansion_intercept
        blanket_expansion_angle = 2. * \
            math.atan(self.blanket_entrance_radius/adjacent)

        self.z_expansion = self.z_entrance + self.dz_expansion
        self.z_active = self.z_expansion + self.dz_active

        # converging region

        self.dz_converging = math.tan(
            self.converging_angle*math.pi/180.)*(self.blanket_active_radius - self.blanket_defueling_radius)
        m = -(self.dz_converging) / \
            (self.seed_active_radius - self.seed_defueling_radius)
        self.seed_converging_intercept = self.z_active - \
            m*self.seed_active_radius
        adjacent = self.z_active - self.seed_converging_intercept
        seed_converging_angle = 2.*math.atan(self.seed_active_radius/adjacent)

        m = -math.tan(self.converging_angle*math.pi/180.)
        self.blanket_converging_intercept = self.z_active - \
            m*self.blanket_active_radius
        adjacent = self.z_active - self.blanket_converging_intercept
        blanket_converging_angle = 2 * \
            math.atan(self.blanket_active_radius/adjacent)

        self.z_converging = self.z_active + self.dz_converging
        self.z_defueling = self.z_converging + self.dz_defueling

        # generate dimensions for shutdown and control rods

        self.inner_shutdown_radius = 2.*math.pi*self.inner_reflector_radius / \
            (4.*(self.inner_shutdown_rods + math.pi))

        # check to see that this is feasible

        if self.inner_reflector_radius < 4.*self.inner_shutdown_radius:
            #print('too few inner shutdown rods - recommend at least 8 shutdown rods')
            poop

        self.outer_reflector_radius = self.blanket_active_radius + \
            self.dr_outer_reflector
        self.outer_shutdown_radius = 2.*math.pi*self.blanket_active_radius / \
            (4.*(self.outer_shutdown_rods/self.shutdown_rod_efficiency - math.pi))

        if self.dr_outer_reflector < 4.*self.outer_shutdown_radius:
            #print('outer reflector is too thin - recommend at least %1.5E cm thickness' % (4.*self.outer_shutdown_radius))
            #print('shutdown rod efficiency should be reduced to %1.5E' % ((2*math.pi*self.blanket_active_radius/self.dr_outer_reflector+math.pi)/self.outer_shutdown_rods))
            #print('number of shutdown rods should be increased to %1.5f' % ((2*math.pi*self.blanket_active_radius/self.dr_outer_reflector+math.pi)/self.shutdown_rod_efficiency))
            poop

        self.dr_shield = self.absorber.pitch
        if self.dr_shield < self.dr_outer_reflector and self.outer_shutdown_rods == 0. or self.dr_shield < (self.dr_outer_reflector - 3.*self.outer_shutdown_radius):
            self.outer_reflector_radius = self.blanket_active_radius + \
                self.dr_outer_reflector - self.dr_shield
            self.shield_radius = self.outer_reflector_radius + self.dr_shield
            # self.no_shield()
        elif self.dr_shield > self.dr_outer_reflector:
            #print('are you trying to simulate the PB-FHR with no shield? if so, have you considered? seting the radius of the absorber pins to something arbitrarily small and keeping the same absorber packing fraction?')
            #print('are you trying to simulate the PB-FHR with no shield? if so, have you considered? using an arbitrary shield parameters (ANNULAR_Graphite.absorber.pitch = 1.) and after you create the input deck changing cell 905 to a solid graphite region?')
            poop
        else:
            #print('shield will cut into shutdown rods')
            #print("have you considered? elimating the outer shutdown rods ('ANNULAR_Graphite.') to 0")
            #print("have you considered? reducing the radius of the outer shutdown rods ('ANNULAR_Graphite.outer_shutdown_radius') to %1.5E" % ((self.dr_outer_reflector - self.dr_shield)/3.))
            #print("have you considered? increasing the graphite reflector radius ('ANNUALR_Graphite.dr_outer_refelctor') to %1.5E" % (3.*self.outer_shutdown_radius + self.dr_shield))
            poop

        self.core_barrel_radius = self.shield_radius + self.dr_core_barrel

        superficial_area = math.pi*(math.pow(self.blanket_active_radius,
                                             2.)-math.pow(self.inner_reflector_radius,
                                                          2.))/math.cos(self.superficial_angle*math.pi/180.)*(1-self.seed.Pebble_PF)
        area = math.pi*math.pow(self.core_barrel_radius, 2.) + superficial_area
        self.down_comer_radius = (
            1. + self.CTE*(self.outlet_temperature - self.temperature_rise/2. - 298.15))*self.shield_radius

        self.pressure_vessel_radius = self.down_comer_radius + \
            self.dr_pressure_vessel

        # generate cells for inner reflector
        theta = 2.*math.pi/self.inner_shutdown_rods
        radius = self.inner_reflector_radius - 2.*self.inner_shutdown_radius
        _Inner_Reflector_Cells_ = '900 1 -1.74 -900 -805 800'
        _Inner_Reflector_Surfaces_ = 'c\nc Inner Reflector Surfaces\nc'
        cell_number = 9000
        _Inner_Shutdown_Channels_ = ''

        for i in range(self.inner_shutdown_rods):
            _Inner_Reflector_Cells_ += '\n          %s' % (cell_number + i)
            angle = i*theta
            x = radius*math.cos(angle)
            y = radius*math.sin(angle)
            _Inner_Reflector_Surfaces_ += '\n%s c/z %1.5E %1.5E %1.5E' % (
                (cell_number + i), x, y, self.inner_shutdown_radius)
            _Inner_Shutdown_Channels_ += '\n%s 24 -_coolant_density_inlet_ -%d -805 800 imp:n=1 tmp=_inlet_coolant_temperature_' % (
                (cell_number + i), (cell_number + i))

        _Inner_Reflector_Cells_ += '\n          imp:n=1 tmp=_inlet_coolant_temperature_' + \
            _Inner_Shutdown_Channels_

        # generate cells for outer reflector
        _Outer_Reflector_Cells_ = '901 1 -1.74 -903 932 -805 800'
        _Outer_Reflector_Surfaces_ = 'c\nc Outer Reflector Surfaces\nc'
        cell_number = 9000 + self.inner_shutdown_rods
        _Outer_Shutdown_Channels_ = ''
        for i in range(self.outer_shutdown_rods):
            theta = 2.*math.pi/self.outer_shutdown_rods
            radius = self.blanket_active_radius + 2*self.outer_shutdown_radius
            _Outer_Reflector_Cells_ += '\n          %s' % (cell_number + i)
            angle = i*theta
            x = radius*math.cos(angle)
            y = radius*math.sin(angle)
            _Outer_Reflector_Surfaces_ += '\n%s c/z %1.5E %1.5E %1.5E' % (
                (cell_number + i), x, y, self.outer_shutdown_radius)
            _Outer_Shutdown_Channels_ += '\n%s 24 -_coolant_density_inlet_ -%d -805 800 imp:n=1 tmp=_inlet_coolant_temperature_' % (
                (cell_number + i), (cell_number + i))

        _Outer_Reflector_Cells_ += '\n          imp:n=1 tmp=_outlet_coolant_temperature_' + \
            _Outer_Shutdown_Channels_

        # generate surfaces for plenum
        dz_plenum = self.plenum_aspect_ratio*self.down_comer_radius
        _plenum1_ = '%1.5E %1.5E %1.5E' % (
            (1. / math.pow(self.pressure_vessel_radius, 2.)),
            (1. / math.pow(self.pressure_vessel_radius, 2.)),
            (1. / math.pow((dz_plenum + self.dr_pressure_vessel), 2.)))
        _plenum2_ = '%1.5E %1.5E %1.5E' % (
            (1. / math.pow(self.down_comer_radius, 2.)),
            (1. / math.pow(self.down_comer_radius, 2.)),
            (1. / math.pow((dz_plenum), 2.)))

        # entrance

        self.dv_seed_entrance = self.dz_active*math.pi * \
            (math.pow(self.seed_entrance_radius, 2.) - math.pow(self.inner_reflector_radius, 2.))
        self.dv_blanket_entrance = self.dz_active*math.pi * \
            (math.pow(self.blanket_entrance_radius, 2.) - math.pow(self.seed_entrance_radius, 2.))

        # seed expansion region volume

        Z = self.z_expansion - self.seed_expansion_intercept
        z = self.z_entrance - self.seed_expansion_intercept

        R = self.seed_active_radius
        r = self.seed_entrance_radius

        dv = math.pi*(1./3.)*math.pow(R, 2.)*Z
        dv += -math.pi*(1./3.)*math.pow(r, 2.)*z
        dv_inner = -dv
        dv += -math.pi * \
            math.pow(self.inner_reflector_radius, 2.)*(self.z_expansion - self.z_entrance)

        self.dv_seed_expansion = dv

        # blanket expansion region volume
        Z = self.z_expansion - self.blanket_expansion_intercept
        z = self.z_entrance - self.blanket_expansion_intercept

        R = self.blanket_active_radius
        r = self.blanket_entrance_radius

        dv = math.pi*(1./3.)*math.pow(R, 2.)*Z
        dv += -math.pi*(1./3.)*math.pow(r, 2.)*z
        dv += dv_inner

        self.dv_blanket_expansion = dv

        # active volume

        self.dv_seed_active = self.dz_active*math.pi * \
            (math.pow(self.seed_active_radius, 2.) - math.pow(self.inner_reflector_radius, 2.))
        self.dv_blanket_active = self.dz_active*math.pi * \
            (math.pow(self.blanket_active_radius, 2.) - math.pow(self.seed_active_radius, 2.))

        # seed converging region volume

        Z = -self.z_active + self.seed_converging_intercept
        z = -self.z_converging + self.seed_converging_intercept

        R = self.seed_active_radius
        r = self.seed_defueling_radius

        dv = math.pi*(1./3.)*math.pow(R, 2.)*Z
        dv += -math.pi*(1./3.)*math.pow(r, 2.)*z
        dv_inner = -dv
        dv += -math.pi * \
            math.pow(self.inner_reflector_radius, 2.)*(self.z_converging - self.z_active)

        self.dv_seed_converging = dv

        # blanket converging region volume

        Z = -self.z_active + self.blanket_converging_intercept
        z = -self.z_converging + self.blanket_converging_intercept

        R = self.blanket_active_radius
        r = self.blanket_defueling_radius

        dv = math.pi*(1./3.)*math.pow(R, 2.)*Z
        dv += -math.pi*(1./3.)*math.pow(r, 2.)*z
        dv += dv_inner

        self.dv_blanket_converging = dv

        self.dv_seed_defueling = self.dz_defueling*math.pi * \
            (math.pow(self.seed_defueling_radius, 2.) - math.pow(self.inner_reflector_radius, 2.))
        self.dv_blanket_defueling = self.dz_defueling*math.pi * \
            (math.pow(self.blanket_defueling_radius, 2.) - math.pow(self.seed_defueling_radius, 2.))

        active_seed_volume = self.dv_seed_expansion + \
            self.dv_seed_active + self.dv_seed_converging
        active_blanket_volume = self.dv_blanket_expansion + \
            self.dv_blanket_active + self.dv_blanket_converging

        self.HT(superficial_area, active_seed_volume)

        active_seed_volume += self.dv_seed_entrance + self.dv_seed_defueling
        active_seed_volume = active_seed_volume*(self.seed.Pebble_PF)*(math.pow(self.seed.r_active,
                                                                                3.)-math.pow(self.seed.r_core,
                                                                                             3.))/math.pow(self.seed.r_shell,
                                                                                                           3.)*(self.seed.TRISO.dv_kernel/self.seed.TRISO.V)/4.

        temperatures = {}
        temperatures['seed'] = self.seed.T_kernel
        temperatures['seed matrix'] = self.seed.T_active
        temperatures['seed core'] = self.seed.T_core
        temperatures['seed shell'] = self.seed.T_shell
        temperatures['coolant'] = self.outlet_temperature - \
            self.temperature_rise/2.
        temperatures['inlet'] = self.outlet_temperature - self.temperature_rise
        temperatures['outlet'] = self.outlet_temperature
        XS = {}

        for index in list(temperatures.keys()):
            if temperatures[index] < (0.5*(293.15 + 600)):
                XS[index] = '70c'
            elif temperatures[index] < 750.:
                XS[index] = '71c'
            elif temperatures[index] < 1050.:
                XS[index] = '72c'
            elif temperatures[index] < (0.5*(1200. + 2500)):
                XS[index] = '73c'
            else:
                XS[index] = '74c'

        _ksrc_card_ = 'ksrc'
        theta = 2*math.pi/self.inner_shutdown_rods
        radius = (self.seed_active_radius + self.inner_reflector_radius)/2.
        z = (self.z_expansion + self.z_active)/2.
        for i in range(self.inner_shutdown_rods):
            angle = i*theta
            x = math.cos(angle)*radius
            y = math.sin(angle)*radius
            _ksrc_card_ += '\n          %1.5E %1.5E %1.5E' % (x, y, z)

        #
        # generate tokens
        #

        self.seed.TRISO.homogenize()
        self.seed.TRISO.matrix()

        _seed_kernel_temperature_ = '%1.5E' % (self.seed.T_kernel*8.617e-11)
        _SeedVolume_ = '%1.5E' % (active_seed_volume)
        _SeedHomoDensity_ = '%1.5E' % self.seed.TRISO.mat_hmatrix.mass()
        _seed_matrix_temperature_ = '%1.5E' % (self.seed.T_active*8.617e-11)

        _seed_pebble_core_density_ = '%1.5E' % self.seed.mat_core.mass()
        _seed_pebble_core_temperature_ = '%1.5E' % (8.617e-11*self.seed.T_core)
        _seed_pebble_shell_density_ = '%1.5E' % self.seed.mat_shell.mass()
        _seed_pebble_shell_temperature_ = '%1.5E' % (
            8.617e-11*self.seed.T_shell)

        _blanket_pebble_shell_density_ = '%1.5E' % self.blanket.mat_shell.mass()
        _blanket_pebble_shell_temperature_ = '%1.5E' % (
            8.617e-11*temperatures['outlet'])

        _coolant_density_ = '%1.5E' % (
            2.28 - .000488
            * (self.outlet_temperature - self.temperature_rise / 2. - 273.15))
        _coolant_temperature_ = '%1.5E' % (
            8.617e-11*(self.outlet_temperature - self.temperature_rise/2.))

        _seed_coolant_density_ = '%1.5E' % (
            2.28 - .000488*(self.seed_T_bulk - 273.15))
        _blanket_coolant_density_ = '%1.5E' % (
            2.28 - .000488*(temperatures['outlet'] - 273.15))
        _seed_coolant_temperature_ = '%1.5E' % (8.617e-11*self.seed_T_bulk)
        _blanket_coolant_temperature_ = '%1.5E' % (
            8.617e-11*temperatures['outlet'])
        _coolant_density_inlet_ = '%1.5E' % (
            2.28 - .000488
            * (self.outlet_temperature - self.temperature_rise - 273.15))
        _inlet_coolant_temperature_ = '%1.5E' % (
            8.617e-11*(self.outlet_temperature - self.temperature_rise))
        _outlet_coolant_temperature_ = '%1.5E' % (
            8.617e-11*(self.outlet_temperature))

        _seed_kernel_radius_ = '%1.5E' % self.seed.TRISO.r_kernel
        _seed_TRISO_hpitch_ = '%1.5E' % (
            math.pow((self.seed.TRISO.V), (1./3.))/2.)

        _seed_pebble_core_ = '%1.5E' % self.seed.r_core
        _seed_pebble_active_ = '%1.5E' % self.seed.r_active
        _seed_pebble_shell_ = '%1.5E' % self.seed.r_shell
        _seed_pebble_hpitch_ = '%1.5E' % self.seed.h_pitch

        _blanket_pebble_shell_ = '%1.5E' % self.blanket.r_shell
        _blanket_pebble_hpitch_ = '%1.5E' % self.blanket.h_pitch

        _inner_reflector_radius_ = '%1.5E' % self.inner_reflector_radius
        _outer_reflector_radius_ = '%1.5E' % self.outer_reflector_radius
        _shield_radius_ = '%1.5E' % self.shield_radius
        _core_barrel_radius_ = '%1.5E' % self.core_barrel_radius
        _down_comer_radius_ = '%1.5E' % self.down_comer_radius
        _pressure_vessel_radius_ = '%1.5E' % self.pressure_vessel_radius

        _absorber_radius_ = '%1.5E' % self.absorber.absorber_radius
        _absorber_hpitch_ = '%1.5E' % (self.absorber.pitch/2.)
        _absorber_pitch_ = '%1.5E' % (self.absorber.pitch)
        _absorber_density_ = '%1.5E' % self.absorber.absorber_density
        _absorber_matrix_density_ = '%1.5E' % self.absorber.mat_matrix.mass()
        _absorber_5010_ = '%1.5E' % (
            self.absorber.mat_absorber.comp['50100']*.60221415)
        _absorber_5011_ = '%1.5E' % (
            self.absorber.mat_absorber.comp['50110']*.60221415)
        _absorber_6012_ = '%1.5E' % (
            self.absorber.mat_absorber.comp['60120']*.60221415)

        _seed_entrance_radius_ = '%1.5E' % self.seed_entrance_radius
        _blanket_entrance_radius_ = '%1.5E' % self.blanket_entrance_radius

        _seed_expansion_tan_ = '%1.5E' % math.pow(
            math.tan(
                seed_expansion_angle/2.),
            2.)
        _seed_expansion_b_ = '%1.5E' % self.seed_expansion_intercept
        _blanket_expansion_tan_ = '%1.5E' % math.pow(
            math.tan(
                blanket_expansion_angle/2.),
            2.)
        _blanket_expansion_b_ = '%1.5E' % self.blanket_expansion_intercept

        _seed_radius_active_ = '%1.5E' % self.seed_active_radius
        _blanket_radius_active_ = '%1.5E' % self.blanket_active_radius

        _seed_converging_tan_ = '%1.5E' % math.pow(
            math.tan(
                seed_converging_angle/2.),
            2.)
        _seed_converging_b_ = '%1.5E' % self.seed_converging_intercept
        _blanket_converging_tan_ = '%1.5E' % math.pow(
            math.tan(
                blanket_converging_angle/2.),
            2.)
        _blanket_converging_b_ = '%1.5E' % self.blanket_converging_intercept

        _seed_radius_defueling_ = '%1.5E' % self.seed_defueling_radius
        _blanket_radius_defueling_ = '%1.5E' % self.blanket_defueling_radius

        _free_surface_Z_ = '%1.5E' % self.z_free_surface
        _entrance_Z_ = '%1.5E' % self.z_entrance
        _expansion_Z_ = '%1.5E' % self.z_expansion
        _active_Z_ = '%1.5E' % self.z_active
        _converging_Z_ = '%1.5E' % self.z_converging
        _defueling_Z_ = '%1.5E' % self.z_defueling

        _seed_carbon_concentration_ = '%1.5E' % (
            self.seed.TRISO.mat_kernel.comp['60120']*.60221415)
        _seed_oxygen_concentration_ = '%1.5E' % (
            self.seed.TRISO.mat_kernel.comp['80160']*.60221415)
        _seed_XS_ = XS['seed']

        _seed_matrix_XS_ = XS['seed matrix']
        _seed_matrix_carbon_concentration_ = '%1.5E' % (
            self.
            seed.
            TRISO.mat_hmatrix.comp
            ['60120'] * .60221415)
        _seed_matrix_silicon_concentration_ = '%1.5E' % (
            self.
            seed.
            TRISO.mat_hmatrix.comp
            ['140280'] * .60221415)

        _coolant_XS_ = XS['coolant']
        _inlet_XS_ = XS['inlet']
        _outlet_XS_ = XS['outlet']
        _seed_core_XS_ = XS['seed core']
        _seed_shell_XS_ = XS['seed shell']
        _blanket_shell_XS_ = XS['outlet']

        _poison_rod_volume_ = '%1.5E' % (
            (math.pi
             *
             (math.pow(self.shield_radius, 2.) -
              math.pow(self.outer_reflector_radius,
                       2.))) * (self.z_defueling) *
            self.absorber.absorber_PF)

        _reflector_surface_ = '%1.5E' % (
            math.pi*2*self.blanket_active_radius*self.dz_active)
        skeleton = 'Full Core PB-AHTR\nc ------------------------------------------------------------------------------\nc Cells\nc ------------------------------------------------------------------------------\nc ------------------------------------------------------------------------------\nc TRISO Particles\nc ------------------------------------------------------------------------------\nc\nc Seed Fuel Particles\nc\nc Depletion Zone 1\n 10  10  4.515775e-02 -1 u=101 tmp=_seed_kernel_temperature_ imp:n=1 $ Kernel 1\n      vol=_SeedVolume_\n11 11 -_SeedHomoDensity_  1 u=101 tmp=_seed_matrix_temperature_ imp:n=1 $ Coatings and Matrix \n12 0 -11 12 -13 14 -15 16  lat=1            imp:n=1 fill=101 u=91 $ TRISO lattice\nc Depletion Zone 2\n 20  20  4.592928e-02 -1 u=102 tmp=_seed_kernel_temperature_ imp:n=1 $ Kernel 2\n      vol=_SeedVolume_\n21 11 -_SeedHomoDensity_  1 u=102 tmp=_seed_matrix_temperature_ imp:n=1 $ Coatings and Matrix \n22 0 -11 12 -13 14 -15 16  lat=1            imp:n=1 fill=102 u=92 $ TRISO lattice\nc Depletion Zone 3\n 30  30  4.676701e-02 -1 u=103 tmp=_seed_kernel_temperature_ imp:n=1 $ Kernel 3\n      vol=_SeedVolume_\n31 11 -_SeedHomoDensity_  1 u=103 tmp=_seed_matrix_temperature_ imp:n=1 $ Coatings and Matrix \n32 0 -11 12 -13 14 -15 16  lat=1            imp:n=1 fill=103 u=93 $ TRISO lattice\nc Depletion Zone 4\n 40  40  4.750444e-02 -1 u=104 tmp=_seed_kernel_temperature_ imp:n=1 $ Kernel 4\n      vol=_SeedVolume_\n41 11 -_SeedHomoDensity_  1 u=104 tmp=_seed_matrix_temperature_ imp:n=1 $ Coatings and Matrix \n42 0 -11 12 -13 14 -15 16  lat=1            imp:n=1 fill=104 u=94 $ TRISO \nc ------------------------------------------------------------------------------\nc Pebbles\nc ------------------------------------------------------------------------------\nc\nc Seed Pebbles\nc\n110 3 -_seed_pebble_core_density_      -110             u=11 imp:n=1 $ Depletion Region 1\n          tmp=_seed_pebble_core_temperature_\n111 0                     -111 110 fill=91 u=11 imp:n=1\n112 4 -_seed_pebble_shell_density_     -112 111         u=11 imp:n=1\n          tmp=_seed_pebble_shell_temperature_\n113 3 -_seed_pebble_core_density_      -113             u=11 imp:n=1 $ Depletion Region 1\n          tmp=_seed_pebble_core_temperature_\n114 0                     -114 113 fill=91 u=11 imp:n=1\n115 4 -_seed_pebble_shell_density_     -115 114         u=11 imp:n=1\n          tmp=_seed_pebble_shell_temperature_\n116 3 -_seed_pebble_core_density_      -116             u=11 imp:n=1 $ Depletion Region 1\n          tmp=_seed_pebble_core_temperature_\n117 0                     -117 116 fill=91 u=11 imp:n=1\n118 4 -_seed_pebble_shell_density_     -118 117         u=11 imp:n=1\n          tmp=_seed_pebble_shell_temperature_\n119 3 -_seed_pebble_core_density_      -119             u=11 imp:n=1 $ Depletion Region 1\n          tmp=_seed_pebble_core_temperature_\n120 0                     -120 119 fill=91 u=11 imp:n=1\n121 4 -_seed_pebble_shell_density_     -121 120         u=11 imp:n=1\n          tmp=_seed_pebble_shell_temperature_\nc\n123 3 -_seed_pebble_core_density_      -130             u=11 imp:n=1 $ Depletion Region 1\n          tmp=_seed_pebble_core_temperature_\n124 0                     -131 130 fill=91 u=11 imp:n=1\n125 4 -_seed_pebble_shell_density_     -132 131         u=11 imp:n=1\n          tmp=_seed_pebble_shell_temperature_\n126 3 -_seed_pebble_core_density_      -133             u=11 imp:n=1 $ Depletion Region 1\n          tmp=_seed_pebble_core_temperature_\n127 0                     -134 133 fill=91 u=11 imp:n=1\n128 4 -_seed_pebble_shell_density_     -135 134         u=11 imp:n=1\n          tmp=_seed_pebble_shell_temperature_\n129 3 -_seed_pebble_core_density_      -136             u=11 imp:n=1 $ Depletion Region 1\n          tmp=_seed_pebble_core_temperature_\n130 0                     -137 136 fill=91 u=11 imp:n=1\n131 4 -_seed_pebble_shell_density_     -138 137         u=11 imp:n=1\n          tmp=_seed_pebble_shell_temperature_\n132 3 -_seed_pebble_core_density_      -139             u=11 imp:n=1 $ Depletion Region 1\n          tmp=_seed_pebble_core_temperature_\n133 0                     -140 139 fill=91 u=11 imp:n=1\n134 4 -_seed_pebble_shell_density_     -141 140         u=11 imp:n=1\n          tmp=_seed_pebble_shell_temperature_\nc\n135 3 -_seed_pebble_core_density_      -150             u=11 imp:n=1 $ Depletion Region 2\n          tmp=_seed_pebble_core_temperature_\n136 0                     -151 150 fill=92 u=11 imp:n=1\n137 4 -_seed_pebble_shell_density_     -152 151         u=11 imp:n=1\n          tmp=_seed_pebble_shell_temperature_\nc\n138 3 -_seed_pebble_core_density_      -153             u=11 imp:n=1 $ Depletion Region 2\n          tmp=_seed_pebble_core_temperature_\n139 0                     -154 153 fill=92 u=11 imp:n=1\n140 4 -_seed_pebble_shell_density_     -155 154         u=11 imp:n=1\n          tmp=_seed_pebble_shell_temperature_\nc\n141 3 -_seed_pebble_core_density_      -160             u=11 imp:n=1 $ Depletion Region 3\n          tmp=_seed_pebble_core_temperature_\n142 0                     -161 160 fill=93 u=11 imp:n=1\n143 4 -_seed_pebble_shell_density_     -162 161         u=11 imp:n=1\n          tmp=_seed_pebble_shell_temperature_\nc\n144 3 -_seed_pebble_core_density_      -163             u=11 imp:n=1 $ Depletion Region 3\n          tmp=_seed_pebble_core_temperature_\n145 0                     -164 163 fill=93 u=11 imp:n=1\n146 4 -_seed_pebble_shell_density_     -165 164         u=11 imp:n=1\n          tmp=_seed_pebble_shell_temperature_\nc\n147 3 -_seed_pebble_core_density_      -170             u=11 imp:n=1 $ Depletion Region 4\n          tmp=_seed_pebble_core_temperature_\n148 0                     -171 170 fill=94 u=11 imp:n=1\n149 4 -_seed_pebble_shell_density_     -172 171         u=11 imp:n=1\n          tmp=_seed_pebble_shell_temperature_\nc\n150 3 -_seed_pebble_core_density_      -173             u=11 imp:n=1 $ Depletion Region 4\n          tmp=_seed_pebble_core_temperature_\n151 0                     -174 173 fill=94 u=11 imp:n=1\n152 4 -_seed_pebble_shell_density_     -175 174         u=11 imp:n=1\n          tmp=_seed_pebble_shell_temperature_\nc\n153 24 -_seed_coolant_density_ 112 115 118 121 \n                   132 135 138 141\n                   152 155 162 165 172 175 u=11 imp:n=1 $ Coolant\n          tmp=_seed_coolant_temperature_\n154 0              -300 301 -302 303 -304 305 fill=11 imp:n=1 u=1 lat=1 $ FCC Unit Cell\nc\nc Blanket Pebbles\nc\n212 6 -_blanket_pebble_shell_density_     -212          u=12 imp:n=1\n          tmp=_blanket_pebble_shell_temperature_\n215 6 -_blanket_pebble_shell_density_     -215          u=12 imp:n=1\n          tmp=_blanket_pebble_shell_temperature_\n218 6 -_blanket_pebble_shell_density_     -218          u=12 imp:n=1\n          tmp=_blanket_pebble_shell_temperature_\n221 6 -_blanket_pebble_shell_density_     -221          u=12 imp:n=1\n          tmp=_blanket_pebble_shell_temperature_\nc\n225 6 -_blanket_pebble_shell_density_     -232          u=12 imp:n=1\n          tmp=_blanket_pebble_shell_temperature_\n228 6 -_blanket_pebble_shell_density_     -235          u=12 imp:n=1\n          tmp=_blanket_pebble_shell_temperature_\n231 6 -_blanket_pebble_shell_density_     -238          u=12 imp:n=1\n          tmp=_blanket_pebble_shell_temperature_\n234 6 -_blanket_pebble_shell_density_     -241          u=12 imp:n=1\n          tmp=_blanket_pebble_shell_temperature_\nc\n237 6 -_blanket_pebble_shell_density_     -252          u=12 imp:n=1\n          tmp=_blanket_pebble_shell_temperature_\nc\n240 6 -_blanket_pebble_shell_density_     -255          u=12 imp:n=1\n          tmp=_blanket_pebble_shell_temperature_\nc\n243 6 -_blanket_pebble_shell_density_     -262          u=12 imp:n=1\n          tmp=_blanket_pebble_shell_temperature_\nc\n246 6 -_blanket_pebble_shell_density_     -265          u=12 imp:n=1\n          tmp=_blanket_pebble_shell_temperature_\nc\n249 6 -_blanket_pebble_shell_density_     -272          u=12 imp:n=1\n          tmp=_blanket_pebble_shell_temperature_\nc\n252 6 -_blanket_pebble_shell_density_     -275          u=12 imp:n=1\n          tmp=_blanket_pebble_shell_temperature_\nc\n253 24 -_blanket_coolant_density_ 212 215 218 221 \n                   232 235 238 241\n                   252 255 262 265 272 275 u=12 imp:n=1 $ Coolant\n          tmp=_blanket_coolant_temperature_\n254 0              -310 311 -312 313 -314 315 fill=12 imp:n=1 u=2 lat=1 $ FCC Unit Cell\nc\nc Radial Shield\nc\n400 400 -_absorber_density_ -400 u=13 imp:n=1 tmp=_outlet_coolant_temperature_ $ b4c poison rod\n          vol=_poison_rod_volume_\n401   2 -_absorber_matrix_density_ 400 u=13 imp:n=1 $ absorber matrix\n402   0 -410 411 -412 413 -414 415 fill=13 lat=2 u=3 imp:n=1 $ hexagonal unit cell \nc\nc Central Graphite Reflector\nc\n_Inner_Reflector_Cells_\n_Outer_Reflector_Cells_\nc\nc lower plenum\nc\n902 24 -_coolant_density_inlet_ -800 -799 imp:n=1 $ lower coolant plenum\n          tmp=_inlet_coolant_temperature_\n903 31 -8.86 -800 -798 799     imp:n=1        $ pressure vessel (plenum)\n          tmp=_inlet_coolant_temperature_\n9040 0       -904 903 -802 800 imp:n=1 fill=3 $ core shield\n9041 0       -904 903 -803 802 imp:n=1 fill=3 $ core shield active\n9042 0       -904 903 -805 803 imp:n=1 fill=3 $ core shield\n905 31 -8.86 -905 904 -805 800 imp:n=1        $ core barrel\n          tmp=_inlet_coolant_temperature_\n906 24 -_coolant_density_inlet_ -906 905 -805 800 imp:n=1 $ down commer\n          tmp=_inlet_coolant_temperature_\n907 31 -8.86 -907 906 -805 800 imp:n=1        $ pressure vessel\n          tmp=_inlet_coolant_temperature_\nc\nc Lower Level\nc\n911 10 -10.5 -911  900 -801 800 imp:n=1 fill=1 $ fuel lower level\n912 20 -9.70 -912  911 -801 800 imp:n=1 fill=2 $ blanket lower level\n913 1  -1.74 -932  912 -801 800 imp:n=1        $ outer radial graphite reflector lower level\n          tmp=_outlet_coolant_temperature_\nc \nc Expansion Region\nc \n921 10 -10.5 -921  900 -802 801 imp:n=1 fill=1 $ fuel expansion region\n922 20 -9.70 -922  921 -802 801 imp:n=1 fill=2 $ blanket expansion region\n923 1  -1.74 -932  922 -802 801 imp:n=1        $ outer radial graphite reflector expansion region\n          tmp=_outlet_coolant_temperature_\nc\nc Active Region\nc\n931 10 -10.5 -931  900 -803 802 imp:n=1 fill=1 $ fuel active region\n932 20 -9.70 -932  931 -803 802 imp:n=1 fill=2 $ blanket active region\nc\nc Converging Region\nc\n941 10 -10.5 -941  900 -804 803 imp:n=1 fill=1 $ fuel converging region\n942 20 -9.70 -942  941 -804 803 imp:n=1 fill=2 $ blanket converging region\n943 2  -1.74 -932  942 -804 803 imp:n=1        $ outer radial graphite reflector converging region\n          tmp=_outlet_coolant_temperature_\nc\nc Upper Level\nc\n951 10 -10.5 -951  900 -805 804 imp:n=1 fill=1 $ fuel upper level\n952 20 -9.70 -952  951 -805 804 imp:n=1 fill=2 $ blanket upper level\n953 2  -1.74 -932  952 -805 804 imp:n=1        $ outer radial graphite reflector upper level\n          tmp=_outlet_coolant_temperature_\nc\n999 0              (800 907:805):(-800 798) imp:n=0        $ void\n\nc ------------------------------------------------------------------------------\nc Surfaces\nc ------------------------------------------------------------------------------\nc TRISO surfaces\nc Seed\n1   so    _seed_kernel_radius_    $  Kernel\nc   Seed Fuel Particles Lattice\n11  px    _seed_TRISO_hpitch_\n12  px   -_seed_TRISO_hpitch_\n13  py    _seed_TRISO_hpitch_\n14  py   -_seed_TRISO_hpitch_\n15  pz    _seed_TRISO_hpitch_\n16  pz   -_seed_TRISO_hpitch_\nc   Blanket Fuel Particles Lattice\nc Seed Pebble Surfaces\n110 s  _seed_pebble_hpitch_  _seed_pebble_hpitch_  _seed_pebble_hpitch_ _seed_pebble_core_ \n111 s  _seed_pebble_hpitch_  _seed_pebble_hpitch_  _seed_pebble_hpitch_ _seed_pebble_active_\n112 s  _seed_pebble_hpitch_  _seed_pebble_hpitch_  _seed_pebble_hpitch_ _seed_pebble_shell_\nc\n113 s -_seed_pebble_hpitch_  _seed_pebble_hpitch_  _seed_pebble_hpitch_ _seed_pebble_core_ \n114 s -_seed_pebble_hpitch_  _seed_pebble_hpitch_  _seed_pebble_hpitch_ _seed_pebble_active_\n115 s -_seed_pebble_hpitch_  _seed_pebble_hpitch_  _seed_pebble_hpitch_ _seed_pebble_shell_\nc\n116 s -_seed_pebble_hpitch_ -_seed_pebble_hpitch_  _seed_pebble_hpitch_ _seed_pebble_core_ \n117 s -_seed_pebble_hpitch_ -_seed_pebble_hpitch_  _seed_pebble_hpitch_ _seed_pebble_active_\n118 s -_seed_pebble_hpitch_ -_seed_pebble_hpitch_  _seed_pebble_hpitch_ _seed_pebble_shell_\nc\n119 s  _seed_pebble_hpitch_ -_seed_pebble_hpitch_  _seed_pebble_hpitch_ _seed_pebble_core_ \n120 s  _seed_pebble_hpitch_ -_seed_pebble_hpitch_  _seed_pebble_hpitch_ _seed_pebble_active_\n121 s  _seed_pebble_hpitch_ -_seed_pebble_hpitch_  _seed_pebble_hpitch_ _seed_pebble_shell_\nc\n130 s  _seed_pebble_hpitch_  _seed_pebble_hpitch_ -_seed_pebble_hpitch_ _seed_pebble_core_ \n131 s  _seed_pebble_hpitch_  _seed_pebble_hpitch_ -_seed_pebble_hpitch_ _seed_pebble_active_ \n132 s  _seed_pebble_hpitch_  _seed_pebble_hpitch_ -_seed_pebble_hpitch_ _seed_pebble_shell_\nc\n133 s -_seed_pebble_hpitch_  _seed_pebble_hpitch_ -_seed_pebble_hpitch_ _seed_pebble_core_\n134 s -_seed_pebble_hpitch_  _seed_pebble_hpitch_ -_seed_pebble_hpitch_ _seed_pebble_active_\n135 s -_seed_pebble_hpitch_  _seed_pebble_hpitch_ -_seed_pebble_hpitch_ _seed_pebble_shell_\nc\n136 s -_seed_pebble_hpitch_ -_seed_pebble_hpitch_ -_seed_pebble_hpitch_ _seed_pebble_core_\n137 s -_seed_pebble_hpitch_ -_seed_pebble_hpitch_ -_seed_pebble_hpitch_ _seed_pebble_active_\n138 s -_seed_pebble_hpitch_ -_seed_pebble_hpitch_ -_seed_pebble_hpitch_ _seed_pebble_shell_\nc\n139 s  _seed_pebble_hpitch_ -_seed_pebble_hpitch_ -_seed_pebble_hpitch_ _seed_pebble_core_\n140 s  _seed_pebble_hpitch_ -_seed_pebble_hpitch_ -_seed_pebble_hpitch_ _seed_pebble_active_\n141 s  _seed_pebble_hpitch_ -_seed_pebble_hpitch_ -_seed_pebble_hpitch_ _seed_pebble_shell_\nc\n150 s  _seed_pebble_hpitch_  0           0          _seed_pebble_core_\n151 s  _seed_pebble_hpitch_  0           0          _seed_pebble_active_\n152 s  _seed_pebble_hpitch_  0           0          _seed_pebble_shell_\nc\n153 s -_seed_pebble_hpitch_  0           0          _seed_pebble_core_\n154 s -_seed_pebble_hpitch_  0           0          _seed_pebble_active_\n155 s -_seed_pebble_hpitch_  0           0          _seed_pebble_shell_\nc\n160 s   0          _seed_pebble_hpitch_  0          _seed_pebble_core_\n161 s   0          _seed_pebble_hpitch_  0          _seed_pebble_active_\n162 s   0          _seed_pebble_hpitch_  0          _seed_pebble_shell_\nc\n163 s   0         -_seed_pebble_hpitch_  0          _seed_pebble_core_\n164 s   0         -_seed_pebble_hpitch_  0          _seed_pebble_active_\n165 s   0         -_seed_pebble_hpitch_  0          _seed_pebble_shell_\nc\n170 s   0           0          _seed_pebble_hpitch_  _seed_pebble_core_\n171 s   0           0          _seed_pebble_hpitch_  _seed_pebble_active_\n172 s   0           0          _seed_pebble_hpitch_  _seed_pebble_shell_\nc\n173 s   0           0         -_seed_pebble_hpitch_  _seed_pebble_core_\n174 s   0           0         -_seed_pebble_hpitch_  _seed_pebble_active_\n175 s   0           0         -_seed_pebble_hpitch_  _seed_pebble_shell_\nc\nc Blanket Pebble Surfaces\n212 s  _blanket_pebble_hpitch_  _blanket_pebble_hpitch_  _blanket_pebble_hpitch_ _blanket_pebble_shell_\nc\n215 s -_blanket_pebble_hpitch_  _blanket_pebble_hpitch_  _blanket_pebble_hpitch_ _blanket_pebble_shell_\nc\n218 s -_blanket_pebble_hpitch_ -_blanket_pebble_hpitch_  _blanket_pebble_hpitch_ _blanket_pebble_shell_\nc\n221 s  _blanket_pebble_hpitch_ -_blanket_pebble_hpitch_  _blanket_pebble_hpitch_ _blanket_pebble_shell_\nc\n232 s  _blanket_pebble_hpitch_  _blanket_pebble_hpitch_ -_blanket_pebble_hpitch_ _blanket_pebble_shell_\nc\n235 s -_blanket_pebble_hpitch_  _blanket_pebble_hpitch_ -_blanket_pebble_hpitch_ _blanket_pebble_shell_\nc\n238 s -_blanket_pebble_hpitch_ -_blanket_pebble_hpitch_ -_blanket_pebble_hpitch_ _blanket_pebble_shell_\nc\n241 s  _blanket_pebble_hpitch_ -_blanket_pebble_hpitch_ -_blanket_pebble_hpitch_ _blanket_pebble_shell_\nc\n252 s  _blanket_pebble_hpitch_  0           0          _blanket_pebble_shell_\nc\n255 s -_blanket_pebble_hpitch_  0           0          _blanket_pebble_shell_\nc\n262 s   0          _blanket_pebble_hpitch_  0          _blanket_pebble_shell_\nc\n265 s   0         -_blanket_pebble_hpitch_  0          _blanket_pebble_shell_\nc\n272 s   0           0          _blanket_pebble_hpitch_  _blanket_pebble_shell_\nc\n275 s   0           0         -_blanket_pebble_hpitch_  _blanket_pebble_shell_\nc\nc Seed Unit Cell Surfaces\n300 px  _seed_pebble_hpitch_\n301 px -_seed_pebble_hpitch_\n302 py  _seed_pebble_hpitch_\n303 py -_seed_pebble_hpitch_\n304 pz  _seed_pebble_hpitch_\n305 pz -_seed_pebble_hpitch_\nc\nc Seed Unit Cell Surfaces\n310 px  _blanket_pebble_hpitch_\n311 px -_blanket_pebble_hpitch_\n312 py  _blanket_pebble_hpitch_\n313 py -_blanket_pebble_hpitch_\n314 pz  _blanket_pebble_hpitch_\n315 pz -_blanket_pebble_hpitch_\nc\nc Shield Absorber Unit Cell Surfaces\n400 cz _absorber_radius_\nc\n410 px  _absorber_hpitch_\n411 px -_absorber_hpitch_\n412 p   1.  1.732050808 0.0  _absorber_pitch_\n413 p   1.  1.732050808 0.0 -_absorber_pitch_\n414 p  -1.  1.732050808 0.0  _absorber_pitch_\n415 p  -1.  1.732050808 0.0 -_absorber_pitch_\nc\n_Inner_Reflector_Surfaces_\n_Outer_Reflector_Surfaces_\nc\n900 cz _inner_reflector_radius_  $ fixed\n903 cz _outer_reflector_radius_\n904 cz _shield_radius_\n905 cz _core_barrel_radius_\n906 cz _down_comer_radius_\n907 cz _pressure_vessel_radius_\nc\nc entrance\nc\n911 cz _seed_entrance_radius_ $ Interface\n912 cz _blanket_entrance_radius_\nc\nc expansion\nc\n921 kz _seed_expansion_b_ _seed_expansion_tan_ 1 $ Interface\n922 kz _blanket_expansion_b_ _blanket_expansion_tan_ 1 $ fixed\nc \nc active\nc\n931 cz _seed_radius_active_ $ Interface\n932 cz _blanket_radius_active_ $ fixed\nc\nc converging\nc\n941 kz _seed_converging_b_ _seed_converging_tan_ -1 $ Interface\n942 kz _blanket_converging_b_ _blanket_converging_tan_ -1 $ fixed\nc\nc defueling chute\nc\n951 cz _seed_radius_defueling_ $ Interface\n952 cz _blanket_radius_defueling_ $ fixed\nc \nc axial subdivision\nc\n798 sq _plenum1_ 0 0 0 -1 0 0 0 $ Lower Surface of Reactor Pressure Vessel\n799 sq _plenum2_ 0 0 0 -1 0 0 0 $ Lower Surface of Coolant Plenum\n800 pz _free_surface_Z_ \n801 pz _entrance_Z_\n802 pz _expansion_Z_\n803 pz _active_Z_\n804 pz _converging_Z_\n805 pz _defueling_Z_\n\nc ------------------------------------------------------------------------------\nc Data\nc ------------------------------------------------------------------------------\nc   Fuel Mix (800 C - Scattering Kernel at 727 C)\nm10        6000._seed_XS_ _seed_carbon_concentration_\n           8016._seed_XS_ _seed_oxygen_concentration_\n          90232._seed_XS_ 3.155364E-12\n          90233._seed_XS_ 6.111943E-18\n          91233._seed_XS_ 2.090064E-13\n          92233._seed_XS_ 1.565088E-11\n          92234._seed_XS_ 2.310757E-08\n          92235._seed_XS_ 6.530410E-03\n          92236._seed_XS_ 2.886212E-04\n          92237._seed_XS_ 6.182108E-07\n          92238._seed_XS_ 3.101251E-02\n          92239._seed_XS_ 5.669246E-08\n          93236._seed_XS_ 6.512238E-13\n          93237._seed_XS_ 7.759563E-06\n          93238._seed_XS_ 2.076991E-08\n          93239._seed_XS_ 8.114427E-06\n          94236._seed_XS_ 1.104723E-12\n          94237._seed_XS_ 2.630942E-13\n          94238._seed_XS_ 6.524353E-07\n          94239._seed_XS_ 2.792373E-04\n          94240._seed_XS_ 4.359378E-05\n          94241._seed_XS_ 2.053468E-05\n          94242._seed_XS_ 1.670386E-06\n          94243._seed_XS_ 4.293403E-10\n          94244._seed_XS_ 1.229758E-11\n          95241._seed_XS_ 2.019546E-07\n          95642._seed_XS_ 7.181079E-09\n          95242._seed_XS_ 4.676686E-10\n          95243._seed_XS_ 9.846851E-08\n          95244._seed_XS_ 4.084775E-12\n          96242._seed_XS_ 2.369615E-08\n          96243._seed_XS_ 1.631013E-10\n          96244._seed_XS_ 8.146734E-09\n          96245._seed_XS_ 2.095970E-10\n          96246._seed_XS_ 3.087975E-12\n          96247._seed_XS_ 1.160300E-14\n          96248._seed_XS_ 2.048319E-16\n          97249._seed_XS_ 1.132133E-18\n          98249._seed_XS_ 7.056902E-20\n          98250._seed_XS_ 2.310151E-19\n          35081._seed_XS_ 2.521353E-06\n          36083._seed_XS_ 5.975651E-06\n          36084._seed_XS_ 1.221480E-05\n          36086._seed_XS_ 2.301418E-05\n          37085._seed_XS_ 1.159947E-05\n          37087._seed_XS_ 2.941082E-05\n          38088._seed_XS_ 4.199715E-05\n          38089._seed_XS_ 1.606480E-05\n          38090._seed_XS_ 6.475389E-05\n          39089._seed_XS_ 3.851261E-05\n          39091._seed_XS_ 2.211768E-05\n          40091._seed_XS_ 4.493853E-05\n          40092._seed_XS_ 6.932220E-05\n          40093._seed_XS_ 7.520800E-05\n          40094._seed_XS_ 7.389051E-05\n          40095._seed_XS_ 2.708780E-05\n          40096._seed_XS_ 7.512218E-05\n          42095._seed_XS_ 3.593013E-05\n          42097._seed_XS_ 7.169469E-05\n          42098._seed_XS_ 7.133124E-05\n          42099._seed_XS_ 1.450855E-06\n          42100._seed_XS_ 7.734324E-05\n          43099._seed_XS_ 7.093751E-05\n          44101._seed_XS_ 6.220977E-05\n          44102._seed_XS_ 5.489038E-05\n          44103._seed_XS_ 1.076203E-05\n          44104._seed_XS_ 2.818823E-05\n          44105._seed_XS_ 2.674051E-08\n          44106._seed_XS_ 8.034671E-06\n          45103._seed_XS_ 2.907716E-05\n          45105._seed_XS_ 1.943475E-07\n          46104._seed_XS_ 3.095244E-06\n          46105._seed_XS_ 1.666650E-05\n          46106._seed_XS_ 3.117808E-06\n          46107._seed_XS_ 6.042788E-06\n          46108._seed_XS_ 3.784125E-06\n          46110._seed_XS_ 1.166105E-06\n          47109._seed_XS_ 2.126308E-06\n          48110._seed_XS_ 2.451390E-07\n          48111._seed_XS_ 6.292656E-07\n          48113._seed_XS_ 2.369867E-08\n          48114._seed_XS_ 5.806548E-07\n          49115._seed_XS_ 1.800116E-07\n          52130._seed_XS_ 1.961243E-05\n          53127._seed_XS_ 2.184106E-06\n          53129._seed_XS_ 9.380429E-06\n          54131._seed_XS_ 3.092417E-05\n          54132._seed_XS_ 5.597567E-05\n          54134._seed_XS_ 9.343075E-05\n          54135._seed_XS_ 8.872110E-08\n          54136._seed_XS_ 1.270242E-04\n          55133._seed_XS_ 7.518780E-05\n          55134._seed_XS_ 2.936893E-06\n          55135._seed_XS_ 3.201097E-05\n          55137._seed_XS_ 7.432967E-05\n          56138._seed_XS_ 8.164906E-05\n          56140._seed_XS_ 6.536468E-06\n          57139._seed_XS_ 7.742905E-05\n          58141._seed_XS_ 1.458123E-05\n          58142._seed_XS_ 7.082646E-05\n          58143._seed_XS_ 6.903447E-07\n          59141._seed_XS_ 5.594538E-05\n          59143._seed_XS_ 6.404719E-06\n          60143._seed_XS_ 5.937792E-05\n          60144._seed_XS_ 2.029238E-05\n          60145._seed_XS_ 4.509401E-05\n          60146._seed_XS_ 3.720168E-05\n          60147._seed_XS_ 2.075930E-06\n          60148._seed_XS_ 2.093547E-05\n          60150._seed_XS_ 8.349153E-06\n          61147._seed_XS_ 1.910310E-05\n          61148._seed_XS_ 1.498052E-07\n          61548._seed_XS_ 1.317339E-07\n          61149._seed_XS_ 2.383294E-07\n          62147._seed_XS_ 1.755089E-06\n          62149._seed_XS_ 7.961982E-07\n          62150._seed_XS_ 1.305224E-05\n          62151._seed_XS_ 2.332160E-06\n          62152._seed_XS_ 5.691457E-06\n          62153._seed_XS_ 7.217423E-08\n          62154._seed_XS_ 1.224509E-06\n          63153._seed_XS_ 3.425777E-06\n          63154._seed_XS_ 3.691345E-07\n          63155._seed_XS_ 2.082947E-07\n          63156._seed_XS_ 1.162370E-07\n          64155._seed_XS_ 2.609943E-09\n          64156._seed_XS_ 7.631348E-07\n          64157._seed_XS_ 1.212848E-08\n          64158._seed_XS_ 3.031187E-07\nc\nc   Fuel Mix (800 C - Scattering Kernel at 727 C)\nm20        6000._seed_XS_ _seed_carbon_concentration_\n           8016._seed_XS_ _seed_oxygen_concentration_\n          90232._seed_XS_ 1.906373E-11\n          90233._seed_XS_ 3.863982E-17\n          91233._seed_XS_ 1.391239E-12\n          92233._seed_XS_ 3.583119E-11\n          92234._seed_XS_ 7.990250E-08\n          92235._seed_XS_ 4.334846E-03\n          92236._seed_XS_ 7.117981E-04\n          92237._seed_XS_ 1.556355E-06\n          92238._seed_XS_ 3.005796E-02\n          92239._seed_XS_ 5.795947E-08\n          93236._seed_XS_ 6.159393E-12\n          93237._seed_XS_ 4.390675E-05\n          93238._seed_XS_ 1.239248E-07\n          93239._seed_XS_ 8.320380E-06\n          94236._seed_XS_ 1.284275E-11\n          94237._seed_XS_ 2.690456E-12\n          94238._seed_XS_ 8.545009E-06\n          94239._seed_XS_ 4.870423E-04\n          94240._seed_XS_ 1.368120E-04\n          94241._seed_XS_ 1.171456E-04\n          94242._seed_XS_ 2.294553E-05\n          94243._seed_XS_ 6.154345E-09\n          94244._seed_XS_ 5.249265E-10\n          95241._seed_XS_ 2.266638E-06\n          95642._seed_XS_ 1.485937E-07\n          95242._seed_XS_ 5.506201E-09\n          95243._seed_XS_ 2.932299E-06\n          95244._seed_XS_ 1.262418E-10\n          96242._seed_XS_ 5.132154E-07\n          96243._seed_XS_ 7.583393E-09\n          96244._seed_XS_ 5.297219E-07\n          96245._seed_XS_ 2.711910E-08\n          96246._seed_XS_ 8.792859E-10\n          96247._seed_XS_ 6.662664E-12\n          96248._seed_XS_ 2.459870E-13\n          97249._seed_XS_ 2.204701E-15\n          98249._seed_XS_ 2.460577E-16\n          98250._seed_XS_ 7.476883E-16\n          35081._seed_XS_ 6.924648E-06\n          36083._seed_XS_ 1.522585E-05\n          36084._seed_XS_ 3.422749E-05\n          36086._seed_XS_ 6.189175E-05\n          37085._seed_XS_ 3.156424E-05\n          37087._seed_XS_ 7.895855E-05\n          38088._seed_XS_ 1.126580E-04\n          38089._seed_XS_ 1.577354E-05\n          38090._seed_XS_ 1.715009E-04\n          39089._seed_XS_ 1.303154E-04\n          39091._seed_XS_ 2.295209E-05\n          40091._seed_XS_ 1.578969E-04\n          40092._seed_XS_ 1.888958E-04\n          40093._seed_XS_ 2.062049E-04\n          40094._seed_XS_ 2.052155E-04\n          40095._seed_XS_ 3.069702E-05\n          40096._seed_XS_ 2.102230E-04\n          42095._seed_XS_ 1.590579E-04\n          42097._seed_XS_ 2.024139E-04\n          42098._seed_XS_ 2.037012E-04\n          42099._seed_XS_ 1.305678E-06\n          42100._seed_XS_ 2.219643E-04\n          43099._seed_XS_ 1.961092E-04\n          44101._seed_XS_ 1.769324E-04\n          44102._seed_XS_ 1.657312E-04\n          44103._seed_XS_ 1.336066E-05\n          44104._seed_XS_ 9.346609E-05\n          44105._seed_XS_ 3.573680E-08\n          44106._seed_XS_ 2.653001E-05\n          45103._seed_XS_ 9.476339E-05\n          45105._seed_XS_ 2.608934E-07\n          46104._seed_XS_ 2.540283E-05\n          46105._seed_XS_ 5.957479E-05\n          46106._seed_XS_ 2.012984E-05\n          46107._seed_XS_ 2.700905E-05\n          46108._seed_XS_ 1.806072E-05\n          46110._seed_XS_ 5.456227E-06\n          47109._seed_XS_ 9.644937E-06\n          48110._seed_XS_ 2.528774E-06\n          48111._seed_XS_ 2.804386E-06\n          48113._seed_XS_ 3.269445E-08\n          48114._seed_XS_ 2.077546E-06\n          49115._seed_XS_ 4.280228E-07\n          52130._seed_XS_ 5.843397E-05\n          53127._seed_XS_ 7.252758E-06\n          53129._seed_XS_ 2.892320E-05\n          54131._seed_XS_ 8.223461E-05\n          54132._seed_XS_ 1.770232E-04\n          54134._seed_XS_ 2.660573E-04\n          54135._seed_XS_ 7.814080E-08\n          54136._seed_XS_ 3.665752E-04\n          55133._seed_XS_ 2.090064E-04\n          55134._seed_XS_ 1.785679E-05\n          55135._seed_XS_ 9.140151E-05\n          55137._seed_XS_ 2.109903E-04\n          56138._seed_XS_ 2.302933E-04\n          56140._seed_XS_ 6.052883E-06\n          57139._seed_XS_ 2.180118E-04\n          58141._seed_XS_ 1.470390E-05\n          58142._seed_XS_ 1.995215E-04\n          58143._seed_XS_ 5.889838E-07\n          59141._seed_XS_ 1.835703E-04\n          59143._seed_XS_ 5.845921E-06\n          60143._seed_XS_ 1.643784E-04\n          60144._seed_XS_ 1.120876E-04\n          60145._seed_XS_ 1.210526E-04\n          60146._seed_XS_ 1.103360E-04\n          60147._seed_XS_ 1.920911E-06\n          60148._seed_XS_ 5.973127E-05\n          60150._seed_XS_ 2.489653E-05\n          61147._seed_XS_ 4.032933E-05\n          61148._seed_XS_ 3.303619E-07\n          61548._seed_XS_ 2.881821E-07\n          61149._seed_XS_ 2.570469E-07\n          62147._seed_XS_ 9.211326E-06\n          62149._seed_XS_ 8.505131E-07\n          62150._seed_XS_ 4.109661E-05\n          62151._seed_XS_ 3.802146E-06\n          62152._seed_XS_ 1.606076E-05\n          62153._seed_XS_ 1.514861E-07\n          62154._seed_XS_ 4.228639E-06\n          63153._seed_XS_ 1.450804E-05\n          63154._seed_XS_ 2.495660E-06\n          63155._seed_XS_ 7.439024E-07\n          63156._seed_XS_ 3.719361E-07\n          64155._seed_XS_ 1.040262E-08\n          64156._seed_XS_ 4.868858E-06\n          64157._seed_XS_ 2.449775E-08\n          64158._seed_XS_ 1.437377E-06\nc\nc   Fuel Mix (800 C - Scattering Kernel at 727 C)\nm30        6000._seed_XS_ _seed_carbon_concentration_\n           8016._seed_XS_ _seed_oxygen_concentration_\n          90232._seed_XS_ 4.455035E-11\n          90233._seed_XS_ 9.062415E-17\n          91233._seed_XS_ 3.077072E-12\n          92233._seed_XS_ 5.268951E-11\n          92234._seed_XS_ 2.266083E-07\n          92235._seed_XS_ 2.843457E-03\n          92236._seed_XS_ 9.561647E-04\n          92237._seed_XS_ 2.038526E-06\n          92238._seed_XS_ 2.910492E-02\n          92239._seed_XS_ 5.701553E-08\n          93236._seed_XS_ 1.567561E-11\n          93237._seed_XS_ 9.172963E-05\n          93238._seed_XS_ 2.605249E-07\n          93239._seed_XS_ 8.187621E-06\n          94236._seed_XS_ 3.938185E-11\n          94237._seed_XS_ 9.392544E-12\n          94238._seed_XS_ 3.024978E-05\n          94239._seed_XS_ 5.150327E-04\n          94240._seed_XS_ 1.754180E-04\n          94241._seed_XS_ 1.984312E-04\n          94242._seed_XS_ 6.957964E-05\n          94243._seed_XS_ 1.827626E-08\n          94244._seed_XS_ 3.509572E-09\n          95241._seed_XS_ 5.616244E-06\n          95642._seed_XS_ 5.067037E-07\n          95242._seed_XS_ 1.367363E-08\n          95243._seed_XS_ 1.418801E-05\n          95244._seed_XS_ 6.109924E-10\n          96242._seed_XS_ 1.805466E-06\n          96243._seed_XS_ 4.222077E-08\n          96244._seed_XS_ 4.290223E-06\n          96245._seed_XS_ 3.249658E-07\n          96246._seed_XS_ 1.762913E-08\n          96247._seed_XS_ 2.058263E-10\n          96248._seed_XS_ 1.212293E-11\n          97249._seed_XS_ 1.391239E-13\n          98249._seed_XS_ 2.178553E-14\n          98250._seed_XS_ 6.057426E-14\n          35081._seed_XS_ 1.031529E-05\n          36083._seed_XS_ 2.095213E-05\n          36084._seed_XS_ 5.260875E-05\n          36086._seed_XS_ 9.115922E-05\n          37085._seed_XS_ 4.699654E-05\n          37087._seed_XS_ 1.160704E-04\n          38088._seed_XS_ 1.655747E-04\n          38089._seed_XS_ 1.159391E-05\n          38090._seed_XS_ 2.481677E-04\n          39089._seed_XS_ 2.025755E-04\n          39091._seed_XS_ 1.719400E-05\n          40091._seed_XS_ 2.493994E-04\n          40092._seed_XS_ 2.812514E-04\n          40093._seed_XS_ 3.080403E-04\n          40094._seed_XS_ 3.110286E-04\n          40095._seed_XS_ 2.459618E-05\n          40096._seed_XS_ 3.203621E-04\n          42095._seed_XS_ 2.680361E-04\n          42097._seed_XS_ 3.099333E-04\n          42098._seed_XS_ 3.161421E-04\n          42099._seed_XS_ 1.082563E-06\n          42100._seed_XS_ 3.453692E-04\n          43099._seed_XS_ 2.888787E-04\n          44101._seed_XS_ 2.718876E-04\n          44102._seed_XS_ 2.706004E-04\n          44103._seed_XS_ 1.240258E-05\n          44104._seed_XS_ 1.617434E-04\n          44105._seed_XS_ 3.711335E-08\n          44106._seed_XS_ 4.110670E-05\n          45103._seed_XS_ 1.401638E-04\n          45105._seed_XS_ 2.723267E-07\n          46104._seed_XS_ 6.793908E-05\n          46105._seed_XS_ 1.067621E-04\n          46106._seed_XS_ 5.204844E-05\n          46107._seed_XS_ 5.468846E-05\n          46108._seed_XS_ 3.777764E-05\n          46110._seed_XS_ 1.146368E-05\n          47109._seed_XS_ 1.876389E-05\n          48110._seed_XS_ 8.346124E-06\n          48111._seed_XS_ 5.806548E-06\n          48113._seed_XS_ 3.574942E-08\n          48114._seed_XS_ 3.802196E-06\n          49115._seed_XS_ 5.344164E-07\n          52130._seed_XS_ 9.335503E-05\n          53127._seed_XS_ 1.205327E-05\n          53129._seed_XS_ 4.662351E-05\n          54131._seed_XS_ 1.123097E-04\n          54132._seed_XS_ 2.972379E-04\n          54134._seed_XS_ 4.122785E-04\n          54135._seed_XS_ 6.505676E-08\n          54136._seed_XS_ 5.738402E-04\n          55133._seed_XS_ 3.066824E-04\n          55134._seed_XS_ 3.958982E-05\n          55135._seed_XS_ 1.422486E-04\n          55137._seed_XS_ 3.245216E-04\n          56138._seed_XS_ 3.541222E-04\n          56140._seed_XS_ 4.935995E-06\n          57139._seed_XS_ 3.342740E-04\n          58141._seed_XS_ 1.198462E-05\n          58142._seed_XS_ 3.058193E-04\n          58143._seed_XS_ 4.693496E-07\n          59141._seed_XS_ 2.909382E-04\n          59143._seed_XS_ 4.669064E-06\n          60143._seed_XS_ 2.306617E-04\n          60144._seed_XS_ 2.388191E-04\n          60145._seed_XS_ 1.758370E-04\n          60146._seed_XS_ 1.789162E-04\n          60147._seed_XS_ 1.581140E-06\n          60148._seed_XS_ 9.258776E-05\n          60150._seed_XS_ 4.001889E-05\n          61147._seed_XS_ 4.376894E-05\n          61148._seed_XS_ 3.576052E-07\n          61548._seed_XS_ 3.121392E-07\n          61149._seed_XS_ 2.347757E-07\n          62147._seed_XS_ 1.739743E-05\n          62149._seed_XS_ 7.773192E-07\n          62150._seed_XS_ 6.531925E-05\n          62151._seed_XS_ 4.342468E-06\n          62152._seed_XS_ 2.190315E-05\n          62153._seed_XS_ 1.912784E-07\n          62154._seed_XS_ 7.532410E-06\n          63153._seed_XS_ 2.682127E-05\n          63154._seed_XS_ 5.687419E-06\n          63155._seed_XS_ 1.526724E-06\n          63156._seed_XS_ 7.255282E-07\n          64155._seed_XS_ 2.367797E-08\n          64156._seed_XS_ 1.433995E-05\n          64157._seed_XS_ 3.965545E-08\n          64158._seed_XS_ 3.297057E-06\nc\nc   Fuel Mix (800 C - Scattering Kernel at 727 C)\nm40        6000._seed_XS_ _seed_carbon_concentration_\n           8016._seed_XS_ _seed_oxygen_concentration_\n          90232._seed_XS_ 7.462244E-11\n          90233._seed_XS_ 1.522282E-16\n          91233._seed_XS_ 4.702178E-12\n          92233._seed_XS_ 7.220452E-11\n          92234._seed_XS_ 5.509229E-07\n          92235._seed_XS_ 1.859226E-03\n          92236._seed_XS_ 1.077465E-03\n          92237._seed_XS_ 2.294149E-06\n          92238._seed_XS_ 2.816602E-02\n          92239._seed_XS_ 5.644512E-08\n          93236._seed_XS_ 2.539828E-11\n          93237._seed_XS_ 1.363022E-04\n          93238._seed_XS_ 3.904920E-07\n          93239._seed_XS_ 8.105846E-06\n          94236._seed_XS_ 7.377945E-11\n          94237._seed_XS_ 2.047763E-11\n          94238._seed_XS_ 6.432987E-05\n          94239._seed_XS_ 5.144269E-04\n          94240._seed_XS_ 1.867807E-04\n          94241._seed_XS_ 2.331554E-04\n          94242._seed_XS_ 1.216634E-04\n          94243._seed_XS_ 3.151578E-08\n          94244._seed_XS_ 1.119513E-08\n          95241._seed_XS_ 8.068492E-06\n          95642._seed_XS_ 8.676254E-07\n          95242._seed_XS_ 1.977851E-08\n          95243._seed_XS_ 3.384738E-05\n          95244._seed_XS_ 1.456357E-09\n          96242._seed_XS_ 3.176615E-06\n          96243._seed_XS_ 1.015174E-07\n          96244._seed_XS_ 1.488108E-05\n          96245._seed_XS_ 1.463777E-06\n          96246._seed_XS_ 1.178523E-07\n          96247._seed_XS_ 1.862709E-09\n          96248._seed_XS_ 1.555547E-10\n          97249._seed_XS_ 2.054023E-12\n          98249._seed_XS_ 4.016578E-13\n          98250._seed_XS_ 1.046269E-12\n          35081._seed_XS_ 1.285537E-05\n          36083._seed_XS_ 2.398943E-05\n          36084._seed_XS_ 6.783813E-05\n          36086._seed_XS_ 1.129306E-04\n          37085._seed_XS_ 5.881761E-05\n          37087._seed_XS_ 1.435105E-04\n          38088._seed_XS_ 2.047158E-04\n          38089._seed_XS_ 8.589430E-06\n          38090._seed_XS_ 3.016649E-04\n          39089._seed_XS_ 2.555527E-04\n          39091._seed_XS_ 1.294371E-05\n          40091._seed_XS_ 3.173587E-04\n          40092._seed_XS_ 3.517800E-04\n          40093._seed_XS_ 3.858833E-04\n          40094._seed_XS_ 3.952774E-04\n          40095._seed_XS_ 1.972853E-05\n          40096._seed_XS_ 4.086087E-04\n          42095._seed_XS_ 3.485796E-04\n          42097._seed_XS_ 3.970340E-04\n          42098._seed_XS_ 4.097899E-04\n          42099._seed_XS_ 9.141666E-07\n          42100._seed_XS_ 4.481688E-04\n          43099._seed_XS_ 3.536931E-04\n          44101._seed_XS_ 3.477518E-04\n          44102._seed_XS_ 3.666106E-04\n          44103._seed_XS_ 1.132587E-05\n          44104._seed_XS_ 2.272393E-04\n          44105._seed_XS_ 3.673980E-08\n          44106._seed_XS_ 4.979962E-05\n          45103._seed_XS_ 1.654889E-04\n          45105._seed_XS_ 2.694444E-07\n          46104._seed_XS_ 1.225316E-04\n          46105._seed_XS_ 1.526169E-04\n          46106._seed_XS_ 9.443023E-05\n          46107._seed_XS_ 8.400641E-05\n          46108._seed_XS_ 5.918105E-05\n          46110._seed_XS_ 1.814250E-05\n          47109._seed_XS_ 2.738209E-05\n          48110._seed_XS_ 1.776946E-05\n          48111._seed_XS_ 9.215365E-06\n          48113._seed_XS_ 3.609721E-08\n          48114._seed_XS_ 5.562736E-06\n          49115._seed_XS_ 5.653598E-07\n          52130._seed_XS_ 1.236018E-04\n          53127._seed_XS_ 1.614809E-05\n          53129._seed_XS_ 6.176051E-05\n          54131._seed_XS_ 1.266154E-04\n          54132._seed_XS_ 4.100625E-04\n          54134._seed_XS_ 5.336593E-04\n          54135._seed_XS_ 5.497619E-08\n          54136._seed_XS_ 7.482941E-04\n          55133._seed_XS_ 3.733646E-04\n          55134._seed_XS_ 6.180594E-05\n          55135._seed_XS_ 1.860589E-04\n          55137._seed_XS_ 4.160240E-04\n          56138._seed_XS_ 4.555386E-04\n          56140._seed_XS_ 4.091791E-06\n          57139._seed_XS_ 4.284821E-04\n          58141._seed_XS_ 9.886225E-06\n          58142._seed_XS_ 3.920568E-04\n          58143._seed_XS_ 3.824558E-07\n          59141._seed_XS_ 3.768527E-04\n          59143._seed_XS_ 3.799117E-06\n          60143._seed_XS_ 2.668852E-04\n          60144._seed_XS_ 3.696797E-04\n          60145._seed_XS_ 2.133526E-04\n          60146._seed_XS_ 2.419790E-04\n          60147._seed_XS_ 1.324355E-06\n          60148._seed_XS_ 1.197503E-04\n          60150._seed_XS_ 5.329021E-05\n          61147._seed_XS_ 4.053478E-05\n          61148._seed_XS_ 3.314371E-07\n          61548._seed_XS_ 2.877177E-07\n          61149._seed_XS_ 2.089610E-07\n          62147._seed_XS_ 2.336248E-05\n          62149._seed_XS_ 6.919095E-07\n          62150._seed_XS_ 8.393574E-05\n          62151._seed_XS_ 4.681381E-06\n          62152._seed_XS_ 2.490309E-05\n          62153._seed_XS_ 2.106520E-07\n          62154._seed_XS_ 1.077667E-05\n          63153._seed_XS_ 3.667671E-05\n          63154._seed_XS_ 8.607098E-06\n          63155._seed_XS_ 2.263206E-06\n          63156._seed_XS_ 1.063078E-06\n          64155._seed_XS_ 3.836976E-08\n          64156._seed_XS_ 2.967937E-05\n          64157._seed_XS_ 5.806043E-08\n          64158._seed_XS_ 5.932239E-06\nc\nc    Single isotopes for depletion\nm101 90232._seed_XS_ 1.\nm102 90233._seed_XS_ 1.\nm103 91233._seed_XS_ 1.\nm104 92233._seed_XS_ 1.\nm105 92234._seed_XS_ 1.\nm106 92235._seed_XS_ 1.\nm107 92236._seed_XS_ 1.\nm108 92237._seed_XS_ 1.\nm109 92238._seed_XS_ 1.\nm110 92239._seed_XS_ 1.\nm111 93236._seed_XS_ 1.\nm112 93237._seed_XS_ 1.\nm113 93238._seed_XS_ 1.\nm114 93239._seed_XS_ 1.\nm115 94236._seed_XS_ 1.\nm116 94237._seed_XS_ 1.\nm117 94238._seed_XS_ 1.\nm118 94239._seed_XS_ 1.\nm119 94240._seed_XS_ 1.\nm120 94241._seed_XS_ 1.\nm121 94242._seed_XS_ 1.\nm122 94243._seed_XS_ 1.\nm123 94244._seed_XS_ 1.\nm124 95241._seed_XS_ 1.\nm125 95642._seed_XS_ 1.\nm126 95242._seed_XS_ 1.\nm127 95243._seed_XS_ 1.\nm128 95244._seed_XS_ 1.\nm129 96242._seed_XS_ 1.\nm130 96243._seed_XS_ 1.\nm131 96244._seed_XS_ 1.\nm132 96245._seed_XS_ 1.\nm133 96246._seed_XS_ 1.\nm134 96247._seed_XS_ 1.\nm135 96248._seed_XS_ 1.\nm136 97249._seed_XS_ 1.\nm137 98249._seed_XS_ 1.\nm138 98250._seed_XS_ 1.\nm139 35081._seed_XS_ 1.\nm140 36083._seed_XS_ 1.\nm141 36084._seed_XS_ 1.\nm142 36086._seed_XS_ 1.\nm143 37085._seed_XS_ 1.\nm144 37087._seed_XS_ 1.\nm145 38088._seed_XS_ 1.\nm146 38089._seed_XS_ 1.\nm147 38090._seed_XS_ 1.\nm148 39089._seed_XS_ 1.\nm149 39091._seed_XS_ 1.\nm150 40091._seed_XS_ 1.\nm151 40092._seed_XS_ 1.\nm152 40093._seed_XS_ 1.\nm153 40094._seed_XS_ 1.\nm154 40095._seed_XS_ 1.\nm155 40096._seed_XS_ 1.\nm156 42095._seed_XS_ 1.\nm157 42097._seed_XS_ 1.\nm158 42098._seed_XS_ 1.\nm159 42099._seed_XS_ 1.\nm160 42100._seed_XS_ 1.\nm161 43099._seed_XS_ 1.\nm162 44101._seed_XS_ 1.\nm163 44102._seed_XS_ 1.\nm164 44103._seed_XS_ 1.\nm165 44104._seed_XS_ 1.\nm166 44105._seed_XS_ 1.\nm167 44106._seed_XS_ 1.\nm168 45103._seed_XS_ 1.\nm169 45105._seed_XS_ 1.\nm170 46104._seed_XS_ 1.\nm171 46105._seed_XS_ 1.\nm172 46106._seed_XS_ 1.\nm173 46107._seed_XS_ 1.\nm174 46108._seed_XS_ 1.\nm175 46110._seed_XS_ 1.\nm176 47109._seed_XS_ 1.\nm177 48110._seed_XS_ 1.\nm178 48111._seed_XS_ 1.\nm179 48113._seed_XS_ 1.\nm180 48114._seed_XS_ 1.\nm181 49115._seed_XS_ 1.\nm182 52130._seed_XS_ 1.\nm183 53127._seed_XS_ 1.\nm184 53129._seed_XS_ 1.\nm185 54131._seed_XS_ 1.\nm186 54132._seed_XS_ 1.\nm187 54134._seed_XS_ 1.\nm188 54135._seed_XS_ 1.\nm189 54136._seed_XS_ 1.\nm190 55133._seed_XS_ 1.\nm191 55134._seed_XS_ 1.\nm192 55135._seed_XS_ 1.\nm193 55137._seed_XS_ 1.\nm194 56138._seed_XS_ 1.\nm195 56140._seed_XS_ 1.\nm196 57139._seed_XS_ 1.\nm197 58141._seed_XS_ 1.\nm198 58142._seed_XS_ 1.\nm199 58143._seed_XS_ 1.\nm200 59141._seed_XS_ 1.\nm201 59143._seed_XS_ 1.\nm202 60143._seed_XS_ 1.\nm203 60144._seed_XS_ 1.\nm204 60145._seed_XS_ 1.\nm205 60146._seed_XS_ 1.\nm206 60147._seed_XS_ 1.\nm207 60148._seed_XS_ 1.\nm208 60150._seed_XS_ 1.\nm209 61147._seed_XS_ 1.\nm210 61148._seed_XS_ 1.\nm211 61548._seed_XS_ 1.\nm212 61149._seed_XS_ 1.\nm213 62147._seed_XS_ 1.\nm214 62149._seed_XS_ 1.\nm215 62150._seed_XS_ 1.\nm216 62151._seed_XS_ 1.\nm217 62152._seed_XS_ 1.\nm218 62153._seed_XS_ 1.\nm219 62154._seed_XS_ 1.\nm220 63153._seed_XS_ 1.\nm221 63154._seed_XS_ 1.\nm222 63155._seed_XS_ 1.\nm223 63156._seed_XS_ 1.\nm224 64155._seed_XS_ 1.\nm225 64156._seed_XS_ 1.\nm226 64157._seed_XS_ 1.\nm227 64158._seed_XS_ 1.\nc\nm1  6000._inlet_XS_  1. $ inlet coolant\nmt1 grph.64t\nm2  6000._outlet_XS_ 1. $ outlet coolant\nmt2 grph.64t\nm3  6000._seed_core_XS_ 1. $ seed pebble core\nmt3 grph.64t\nm4  6000._seed_shell_XS_ 1. $ seed pebble shell\nmt4 grph.64t\nm6  6000._blanket_shell_XS_   1. $ blanket pebble shell\nmt6 grph.64t\nc   Coatings and matrix mix\nm11    6000._seed_matrix_XS_  _seed_matrix_carbon_concentration_\n      14028._seed_matrix_XS_  _seed_matrix_silicon_concentration_\nmt11    grph.64t\nc   Salt in the Core - 2LiF-BeF2 - (655 C)\nm24    3006._coolant_XS_  4.12E-07\n       3007._coolant_XS_  2.38358E-02\n       4009._coolant_XS_  1.19185E-02\n       9019._coolant_XS_  4.76740E-02\nm536   3006._coolant_XS_  1\nm549   4009._coolant_XS_  1\nc Steel for Reactor Pressure Vessel and Core Barrel\nm31       28058._inlet_XS_ 0.0438374141756552\n          28060._inlet_XS_ 0.0163231374502\n          28061._inlet_XS_ 0.00069798731409836\n          28062._inlet_XS_ 0.00218910073380645\n          28064._inlet_XS_ 0.000540385299125\n          24050._inlet_XS_ 0.00032455797836\n          24052._inlet_XS_ 0.00601805533522308\n          24053._inlet_XS_ 0.000669523710271698\n          24054._inlet_XS_ 0.000163572351962963\n          42092._inlet_XS_ 0.00154914850330435\n          42094._inlet_XS_ 0.000945063210638298\n          42095._inlet_XS_ 0.00160940903949474\n          42096._inlet_XS_ 0.001668675123\n          42097._inlet_XS_ 0.000945537190515464\n          42098._inlet_XS_ 0.00236471183191837\n          42100._inlet_XS_ 0.00092485418328\n          26054._inlet_XS_ 0.000286535681481481\n          26056._inlet_XS_ 0.00436938684142857\n          26057._inlet_XS_ 0.000102965635087719\n          26058._inlet_XS_ 1.28787737931034E-05\nm400      5010._outlet_XS_ _absorber_5010_\n          5011._outlet_XS_ _absorber_5011_\n          6000._outlet_XS_ _absorber_6012_\nfc34 tallies for boron shielding\n f34:n 400<9041\nfm34      (1)\n          (1 400 103:104:105:106:107)\n         (-1 400 -2)\nfc44 tallies for radiation damage in core barrel\n f44:n 904\nfm44      (1)\n          (1 1 444)\n e44 .1 20\nfc54 tallies for radiation damage in RPV\n f54:n 907\nfm54      (1)\n          (-1 31 444)\n          (1 31 103:104:105:106:107)\n e54 .1 20\n f64:n (153 253 )\nfm64 (1 536 -2)\nsd64 1\n f74:n (153 253 )\nfm74 (1 549 107)\nsd74 1\nc    begin_mocup_reaction_rate_tallies\nc    time dependent reaction rates\nfc104  Reaction rates\nf104:n   \n       10 20 30 40  \nfm104  (1)\n       (1 101 (16) (17) (-6) (102) )\n       (1 102 (16) (17) (-6) (102) )\n       (1 103 (16) (17) (-6) (102) )\n       (1 104 (16) (17) (-6) (102) )\n       (1 105 (16) (17) (-6) (102) )\n       (1 106 (16) (17) (-6) (102) )\n       (1 107 (16) (17) (-6) (102) )\n       (1 108 (16) (17) (-6) (102) )\n       (1 109 (16) (17) (-6) (102) )\n       (1 110 (16) (17) (-6) (102) )\n       (1 111 (16) (17) (-6) (102) )\n       (1 112 (16) (17) (-6) (102) )\n       (1 113 (16) (17) (-6) (102) )\n       (1 114 (16) (17) (-6) (102) )\n       (1 115 (16) (17) (-6) (102) )\n       (1 116 (16) (17) (-6) (102) )\n       (1 117 (16) (17) (-6) (102) )\n       (1 118 (16) (17) (-6) (102) )\n       (1 119 (16) (17) (-6) (102) )\n       (1 120 (16) (17) (-6) (102) )\n       (1 121 (16) (17) (-6) (102) )\n       (1 122 (16) (17) (-6) (102) )\n       (1 123 (16) (17) (-6) (102) )\n       (1 124 (16) (17) (-6) (102) )\n       (1 125 (16) (17) (-6) (102) )\n       (1 126 (16) (17) (-6) (102) )\n       (1 127 (16) (17) (-6) (102) )\n       (1 128 (16) (17) (-6) (102) )\n       (1 129 (16) (17) (-6) (102) )\n       (1 130 (16) (17) (-6) (102) )\n       (1 131 (16) (17) (-6) (102) )\n       (1 132 (16) (17) (-6) (102) )\n       (1 133 (16) (17) (-6) (102) )\n       (1 134 (16) (17) (-6) (102) )\n       (1 135 (16) (17) (-6) (102) )\n       (1 136 (16) (17) (-6) (102) )\n       (1 137 (16) (17) (-6) (102) )\n       (1 138 (16) (17) (-6) (102) )\n       (1 139 (16) (17)      (102) )\n       (1 140 (16) (17)      (102) )\n       (1 141 (16) (17)      (102) )\n       (1 142 (16) (17)      (102) )\n       (1 143 (16) (17)      (102) )\n       (1 144 (16) (17)      (102) )\n       (1 145 (16) (17)      (102) )\n       (1 146 (16) (17)      (102) )\n       (1 147 (16) (17)      (102) )\n       (1 148 (16) (17)      (102) )\n       (1 149 (16) (17)      (102) )\n       (1 150 (16) (17)      (102) )\n       (1 151 (16) (17)      (102) )\n       (1 152 (16) (17)      (102) )\n       (1 153 (16) (17)      (102) )\n       (1 154 (16) (17)      (102) )\n       (1 155 (16) (17)      (102) )\n       (1 156 (16) (17)      (102) )\n       (1 157 (16) (17)      (102) )\n       (1 158 (16) (17)      (102) )\n       (1 159 (16) (17)      (102) )\n       (1 160 (16) (17)      (102) )\n       (1 161 (16) (17)      (102) )\n       (1 162 (16) (17)      (102) )\n       (1 163 (16) (17)      (102) )\n       (1 164 (16) (17)      (102) )\n       (1 165 (16) (17)      (102) )\n       (1 166 (16) (17)      (102) )\n       (1 167 (16) (17)      (102) )\n       (1 168 (16) (17)      (102) )\n       (1 169 (16) (17)      (102) )\n       (1 170 (16) (17)      (102) )\n       (1 171 (16) (17)      (102) )\n       (1 172 (16) (17)      (102) )\n       (1 173 (16) (17)      (102) )\n       (1 174 (16) (17)      (102) )\n       (1 175 (16) (17)      (102) )\n       (1 176 (16) (17)      (102) )\n       (1 177 (16) (17)      (102) )\n       (1 178 (16) (17)      (102) )\n       (1 179 (16) (17)      (102) )\n       (1 180 (16) (17)      (102) )\n       (1 181 (16) (17)      (102) )\n       (1 182 (16) (17)      (102) )\n       (1 183 (16) (17)      (102) )\n       (1 184 (16) (17)      (102) )\n       (1 185 (16) (17)      (102) )\n       (1 186 (16) (17)      (102) )\n       (1 187 (16) (17)      (102) )\n       (1 188 (16) (17)      (102) )\n       (1 189 (16) (17)      (102) )\n       (1 190 (16) (17)      (102) )\n       (1 191 (16) (17)      (102) )\n       (1 192 (16) (17)      (102) )\n       (1 193 (16) (17)      (102) )\n       (1 194 (16) (17)      (102) )\n       (1 195 (16) (17)      (102) )\n       (1 196 (16) (17)      (102) )\n       (1 197 (16) (17)      (102) )\n       (1 198 (16) (17)      (102) )\n       (1 199 (16) (17)      (102) )\n       (1 200 (16) (17)      (102) )\n       (1 201 (16) (17)      (102) )\n       (1 202 (16) (17)      (102) )\n       (1 203 (16) (17)      (102) )\n       (1 204 (16) (17)      (102) )\n       (1 205 (16) (17)      (102) )\n       (1 206 (16) (17)      (102) )\n       (1 207 (16) (17)      (102) )\n       (1 208 (16) (17)      (102) )\n       (1 209 (16) (17)      (102) )\n       (1 210 (16) (17)      (102) )\n       (1 211 (16) (17)      (102) )\n       (1 212 (16) (17)      (102) )\n       (1 213 (16) (17)      (102) )\n       (1 214 (16) (17)      (102) )\n       (1 215 (16) (17)      (102) )\n       (1 216 (16) (17)      (102) )\n       (1 217 (16) (17)      (102) )\n       (1 218 (16) (17)      (102) )\n       (1 219 (16) (17)      (102) )\n       (1 220 (16) (17)      (102) )\n       (1 221 (16) (17)      (102) )\n       (1 222 (16) (17)      (102) )\n       (1 223 (16) (17)      (102) )\n       (1 224 (16) (17)      (102) )\n       (1 225 (16) (17)      (102) )\n       (1 226 (16) (17)      (102) )\n       (1 227 (16) (17)      (102) )\nc\nfc2 tally of fluence on first wall of outer graphite reflector\nf2:n 932\nfs2 -802 803\n e2 .1 20\nsd2 1 1 _reflector_surface_\nc\nfc404 Neutron Spectrum Tally\nf404:n 10 20 30 40 \ne404 1.0000E-11 1.0000E-10 5.0000E-10 7.5000E-10 1.0000E-09 1.2000E-09\n     1.5000E-09 2.0000E-09 2.5000E-09 3.0000E-09 4.0000E-09 5.0000E-09\n     7.5000E-09 1.0000E-08 2.5300E-08 3.0000E-08 4.0000E-08 5.0000E-08\n     6.0000E-08 7.0000E-08 8.0000E-08 9.0000E-08 1.0000E-07 1.2500E-07\n     1.5000E-07 1.7500E-07 2.0000E-07 2.2500E-07 2.5000E-07 2.7500E-07\n     3.0000E-07 3.2500E-07 3.5000E-07 3.7500E-07 4.0000E-07 4.5000E-07\n     5.0000E-07 5.5000E-07 6.0000E-07 6.2500E-07 6.5000E-07 7.0000E-07\n     7.5000E-07 8.0000E-07 8.5000E-07 9.0000E-07 9.2500E-07 9.5000E-07\n     9.7500E-07 1.0000E-06 1.0100E-06 1.0200E-06 1.0300E-06 1.0400E-06\n     1.0500E-06 1.0600E-06 1.0700E-06 1.0800E-06 1.0900E-06 1.1000E-06\n     1.1100E-06 1.1200E-06 1.1300E-06 1.1400E-06 1.1500E-06 1.1750E-06\n     1.2000E-06 1.2250E-06 1.2500E-06 1.3000E-06 1.3500E-06 1.4000E-06\n     1.4500E-06 1.5000E-06 1.5900E-06 1.6800E-06 1.7700E-06 1.8600E-06\n     1.9400E-06 2.0000E-06 2.1200E-06 2.2100E-06 2.3000E-06 2.3800E-06\n     2.4700E-06 2.5700E-06 2.6700E-06 2.7700E-06 2.8700E-06 2.9700E-06\n     3.0000E-06 3.0500E-06 3.1500E-06 3.5000E-06 3.7300E-06 4.0000E-06\n     4.7500E-06 5.0000E-06 5.4000E-06 6.0000E-06 6.2500E-06 6.5000E-06\n     6.7500E-06 7.0000E-06 7.1500E-06 8.1000E-06 9.1000E-06 1.0000E-05\n     1.1500E-05 1.1900E-05 1.2900E-05 1.3750E-05 1.4400E-05 1.5100E-05\n     1.6000E-05 1.7000E-05 1.8500E-05 1.9000E-05 2.0000E-05 2.1000E-05\n     2.2500E-05 2.5000E-05 2.7500E-05 3.0000E-05 3.1250E-05 3.1750E-05\n     3.3250E-05 3.3750E-05 3.4600E-05 3.5500E-05 3.7000E-05 3.8000E-05\n     3.9100E-05 3.9600E-05 4.1000E-05 4.2400E-05 4.4000E-05 4.5200E-05\n     4.7000E-05 4.8300E-05 4.9200E-05 5.0600E-05 5.2000E-05 5.3400E-05\n     5.9000E-05 6.1000E-05 6.5000E-05 6.7500E-05 7.2000E-05 7.6000E-05\n     8.0000E-05 8.2000E-05 9.0000E-05 1.0000E-04 1.0800E-04 1.1500E-04\n     1.1900E-04 1.2200E-04 1.8600E-04 1.9250E-04 2.0750E-04 2.1000E-04\n     2.4000E-04 2.8500E-04 3.0500E-04 5.5000E-04 6.7000E-04 6.8300E-04\n     9.5000E-04 1.1500E-03 1.5000E-03 1.5500E-03 1.8000E-03 2.2000E-03\n     2.2900E-03 2.5800E-03 3.0000E-03 3.7400E-03 3.9000E-03 6.0000E-03\n     8.0300E-03 9.5000E-03 1.3000E-02 1.7000E-02 2.5000E-02 3.0000E-02\n     4.5000E-02 5.0000E-02 5.2000E-02 6.0000E-02 7.3000E-02 7.5000E-02\n     8.2000E-02 8.5000E-02 1.0000E-01 1.2830E-01 1.5000E-01 2.0000E-01\n     2.7000E-01 3.3000E-01 4.0000E-01 4.2000E-01 4.4000E-01 4.7000E-01\n     4.9952E-01 5.5000E-01 5.7300E-01 6.0000E-01 6.7000E-01 6.7900E-01\n     7.5000E-01 8.2000E-01 8.6110E-01 8.7500E-01 9.0000E-01 9.2000E-01\n     1.0100E+00 1.1000E+00 1.2000E+00 1.2500E+00 1.3170E+00 1.3560E+00\n     1.4000E+00 1.5000E+00 1.8500E+00 2.3540E+00 2.4790E+00 3.0000E+00\n     4.3040E+00 4.8000E+00 6.4340E+00 8.1873E+00 1.0000E+01 1.2840E+01\n     1.3840E+01 1.4550E+01 1.5683E+01 1.7333E+01 2.0000E+01\nsd404 1 3r   \nkcode  10000  1.0 10 110\nprdmp  10000 10000 10000\n#print\nmode n\n_ksrc_card_\n\n\n\n\n\n\n\n\n\n\n'

        # replace values

        skeleton = skeleton.replace(
            '_Outer_Reflector_Cells_',
            _Outer_Reflector_Cells_)
        skeleton = skeleton.replace(
            '_Outer_Reflector_Surfaces_',
            _Outer_Reflector_Surfaces_)
        skeleton = skeleton.replace(
            '_Inner_Reflector_Cells_',
            _Inner_Reflector_Cells_)
        skeleton = skeleton.replace(
            '_Inner_Reflector_Surfaces_',
            _Inner_Reflector_Surfaces_)
        skeleton = skeleton.replace(
            '_seed_kernel_temperature_',
            _seed_kernel_temperature_)
        skeleton = skeleton.replace('_SeedVolume_', _SeedVolume_)
        skeleton = skeleton.replace('_SeedHomoDensity_', _SeedHomoDensity_)
        skeleton = skeleton.replace(
            '_seed_matrix_temperature_',
            _seed_matrix_temperature_)

        skeleton = skeleton.replace(
            '_seed_pebble_core_density_',
            _seed_pebble_core_density_)
        skeleton = skeleton.replace(
            '_seed_pebble_core_temperature_',
            _seed_pebble_core_temperature_)
        skeleton = skeleton.replace(
            '_seed_pebble_shell_density_',
            _seed_pebble_shell_density_)
        skeleton = skeleton.replace(
            '_seed_pebble_shell_temperature_',
            _seed_pebble_shell_temperature_)
        skeleton = skeleton.replace(
            '_blanket_pebble_shell_density_',
            _blanket_pebble_shell_density_)
        skeleton = skeleton.replace(
            '_blanket_pebble_shell_temperature_',
            _blanket_pebble_shell_temperature_)

        skeleton = skeleton.replace(
            '_inlet_coolant_temperature_',
            _inlet_coolant_temperature_)
        skeleton = skeleton.replace(
            '_coolant_density_inlet_',
            _coolant_density_inlet_)
        skeleton = skeleton.replace(
            '_outlet_coolant_temperature_',
            _outlet_coolant_temperature_)
        skeleton = skeleton.replace(
            '_seed_coolant_density_',
            _seed_coolant_density_)
        skeleton = skeleton.replace(
            '_blanket_coolant_density_',
            _blanket_coolant_density_)
        skeleton = skeleton.replace(
            '_seed_coolant_temperature_',
            _seed_coolant_temperature_)
        skeleton = skeleton.replace(
            '_blanket_coolant_temperature_',
            _blanket_coolant_temperature_)
        skeleton = skeleton.replace(
            '_coolant_temperature_',
            _coolant_temperature_)
        skeleton = skeleton.replace('_coolant_density_', _coolant_density_)
        skeleton = skeleton.replace(
            '_seed_kernel_radius_',
            _seed_kernel_radius_)
        skeleton = skeleton.replace('_seed_TRISO_hpitch_', _seed_TRISO_hpitch_)
        skeleton = skeleton.replace('_seed_pebble_core_', _seed_pebble_core_)

        skeleton = skeleton.replace(
            '_seed_pebble_active_',
            _seed_pebble_active_)
        skeleton = skeleton.replace('_seed_pebble_shell_', _seed_pebble_shell_)
        skeleton = skeleton.replace(
            '_seed_pebble_hpitch_',
            _seed_pebble_hpitch_)
        skeleton = skeleton.replace(
            '_blanket_pebble_shell_',
            _blanket_pebble_shell_)
        skeleton = skeleton.replace(
            '_blanket_pebble_hpitch_',
            _blanket_pebble_hpitch_)
        skeleton = skeleton.replace(
            '_inner_reflector_radius_',
            _inner_reflector_radius_)
        skeleton = skeleton.replace(
            '_outer_reflector_radius_',
            _outer_reflector_radius_)
        skeleton = skeleton.replace('_shield_radius_', _shield_radius_)

        skeleton = skeleton.replace(
            '_core_barrel_radius_',
            _core_barrel_radius_)
        skeleton = skeleton.replace('_down_comer_radius_', _down_comer_radius_)
        skeleton = skeleton.replace(
            '_pressure_vessel_radius_',
            _pressure_vessel_radius_)
        skeleton = skeleton.replace('_shield_radius_', _shield_radius_)
        skeleton = skeleton.replace('_absorber_radius_', _absorber_radius_)
        skeleton = skeleton.replace('_absorber_hpitch_', _absorber_hpitch_)
        skeleton = skeleton.replace('_absorber_pitch_', _absorber_pitch_)
        skeleton = skeleton.replace('_absorber_density_', _absorber_density_)
        skeleton = skeleton.replace(
            '_absorber_matrix_density_',
            _absorber_matrix_density_)
        skeleton = skeleton.replace('_absorber_5010_', _absorber_5010_)
        skeleton = skeleton.replace('_absorber_5011_', _absorber_5011_)
        skeleton = skeleton.replace('_absorber_6012_', _absorber_6012_)

        skeleton = skeleton.replace(
            '_seed_entrance_radius_',
            _seed_entrance_radius_)
        skeleton = skeleton.replace(
            '_blanket_entrance_radius_',
            _blanket_entrance_radius_)
        skeleton = skeleton.replace(
            '_seed_expansion_tan_',
            _seed_expansion_tan_)
        skeleton = skeleton.replace('_seed_expansion_b_', _seed_expansion_b_)
        skeleton = skeleton.replace(
            '_blanket_expansion_b_',
            _blanket_expansion_b_)
        skeleton = skeleton.replace(
            '_blanket_expansion_tan_',
            _blanket_expansion_tan_)
        skeleton = skeleton.replace(
            '_seed_radius_active_',
            _seed_radius_active_)
        skeleton = skeleton.replace(
            '_blanket_radius_active_',
            _blanket_radius_active_)
        skeleton = skeleton.replace('_seed_converging_b_', _seed_converging_b_)
        skeleton = skeleton.replace(
            '_seed_converging_tan_',
            _seed_converging_tan_)

        skeleton = skeleton.replace(
            '_blanket_converging_b_',
            _blanket_converging_b_)
        skeleton = skeleton.replace(
            '_blanket_converging_tan_',
            _blanket_converging_tan_)
        skeleton = skeleton.replace(
            '_seed_radius_defueling_',
            _seed_radius_defueling_)
        skeleton = skeleton.replace(
            '_blanket_radius_defueling_',
            _blanket_radius_defueling_)
        skeleton = skeleton.replace('_plenum1_', _plenum1_)
        skeleton = skeleton.replace('_plenum2_', _plenum2_)
        skeleton = skeleton.replace('_free_surface_Z_', _free_surface_Z_)
        skeleton = skeleton.replace('_entrance_Z_', _entrance_Z_)
        skeleton = skeleton.replace('_expansion_Z_', _expansion_Z_)
        skeleton = skeleton.replace('_active_Z_', _active_Z_)

        skeleton = skeleton.replace('_converging_Z_', _converging_Z_)
        skeleton = skeleton.replace('_defueling_Z_', _defueling_Z_)
        skeleton = skeleton.replace(
            '_seed_carbon_concentration_',
            _seed_carbon_concentration_)
        skeleton = skeleton.replace(
            '_seed_oxygen_concentration_',
            _seed_oxygen_concentration_)
        skeleton = skeleton.replace('_seed_XS_', _seed_XS_)
        skeleton = skeleton.replace('_seed_matrix_XS_', _seed_matrix_XS_)
        skeleton = skeleton.replace(
            '_seed_matrix_carbon_concentration_',
            _seed_matrix_carbon_concentration_)
        skeleton = skeleton.replace(
            '_seed_matrix_silicon_concentration_',
            _seed_matrix_silicon_concentration_)

        skeleton = skeleton.replace('_coolant_XS_', _coolant_XS_)
        skeleton = skeleton.replace('_inlet_XS_', _inlet_XS_)
        skeleton = skeleton.replace('_outlet_XS_', _outlet_XS_)
        skeleton = skeleton.replace('_seed_core_XS_', _seed_core_XS_)
        skeleton = skeleton.replace('_seed_shell_XS_', _seed_shell_XS_)
        skeleton = skeleton.replace('_blanket_shell_XS_', _blanket_shell_XS_)
        skeleton = skeleton.replace('_reflector_surface_', _reflector_surface_)

        skeleton = skeleton.replace('_poison_rod_volume_', _poison_rod_volume_)

        self.input_source = skeleton.replace('_ksrc_card_', _ksrc_card_)
        self.input = skeleton.replace('_ksrc_card_', '')

        open('inp.s', 'w').write(self.input_source)

    def HT(self, superficial_area, seed_volume):

        import math

        # convert areas and volumes in terms of meters
        superficial_area = superficial_area*1e-4
        seed_volume = seed_volume*1e-6

        bulk_T = self.outlet_temperature - self.temperature_rise/2.
        seed_T_bulk = bulk_T

        pCp = (2282 - .488*(bulk_T - 273.15))*2368
        # volumetric heat capacity in (j/m3*K)
        p = (2282 - .488*(bulk_T - 273.15))
        # density (kg/m3)
        u = .116*.001*math.exp(3755./bulk_T)
        # viscosity (kg/m*s)

        power = (seed_volume)*self.power_density
        superficial_velocity = power*1.e+06 / \
            (pCp*self.temperature_rise*superficial_area)

        # calculate Re

        Re = p*superficial_velocity*2*self.seed.r_shell/u

        seed_power = power

        seed_power_density = seed_power/seed_volume

        self.seed.HT(seed_T_bulk, Re, seed_power_density)
        self.seed_T_bulk = seed_T_bulk

    def DECF(self):

        # This method generates an MCNP5 input deck that can be used on the DECF
        # cluster

        # generate input deck

        self.generate()

        DECF = self.input

        # remove ENDFVII cross sections

        DECF = DECF.replace('.71c', '')
        DECF = DECF.replace('.72c', '')
        DECF = DECF.replace('.73c', '')

        # replace meta stable isotopes

        DECF = DECF.replace('95642', '95601')
        DECF = DECF.replace('61548', '61601')

        # remove pebble heterogeneity (so you can quickly use mcnp plotter)

        DECF = DECF.replace('fill=1 ', '')
        DECF = DECF.replace('fill=2 ', '')

        # add AW table
        awtab = 'awtab  36085  84.183100  38089  88.143700  38090  89.135400\n       44105 104.007000  44106 104.997000  46107 105.986755\n       50123 121.850000  50126 124.826000  51124 122.842000\n       51125 123.832000  52601 125.815000  52129 127.800000\n       54133 131.764160  57140 138.708000  58141 139.698000\n       58143 141.685000  58144 142.677000  61601 146.646787\n       59142 140.691000  59143 141.683000  62153 151.608000\n       62154 152.600000  63156 154.585000  65160 157.562000\n'

        k = DECF.find('m10 ')
        DECF = DECF[:k] + awtab + DECF[k:]
        k = DECF.find('ksrc')
        DECF = DECF[:k] + 'void\n' + DECF[k:]

        return DECF

    def BEAU(self):

        # This method sets up a BEAU equilibrium depletion analysis

        import os
        self.generate()

        skeleton = "#!/usr/bin/env python\n\nimport BEAU\nimport mocup\n\n# -----------------------------------------------------------------------------\n# DEFINE EQUILIBRIUM DEPLETION ANALYSIS PROBLEM\n# -----------------------------------------------------------------------------\n\n# initialize Equilibrium Depletion\nEq_Cycle = BEAU.equilibrium_cycle()\n\n# DEFINE POWER OF SYSTEM IN MWth\nEq_Cycle.power = _power_\n\t# 8 1/2 pebbles with a volume of 1.17810E+01 cc with a power desnity of 30 MW/m3\n\n# DEFINE LOADING PATTER\n\t# 'Y':'X' EOEC material from X is advanced to cell Y\nEq_Cycle.loading_pattern = {\n\t\t\t\t\t\t\t'20':'10',\n\t\t\t\t\t\t\t'30':'20',\n\t\t\t\t\t\t\t'40':'30',\n\t\t\t\t\t\t\t}\n\n# DEFINE THE BURNUPS OF ALL FUEL PROGRESSIONS\nEq_Cycle.burnup = {'seed': 1.75E+5,}\n\n# DEFINE THE MAKE UP MATERIAL FOR EACH PROGRESSION\n\n# initiate material object\nseed = mocup.material()\n# define its composition vector\n\nseed.comp['922350'] = _molar_mass_92235_\nseed.comp['922380'] = _molar_mass_92238_\n\n# link material object to the progression key\nEq_Cycle.makeup = {'seed':seed,}\n\n# INITIATE A PROGRESSION\nseed_progression = Eq_Cycle.fuel_progression()\n\n# DEFINE THE ORDER IN WHICH FUEL IS ADVANCED\nseed_progression.cells = ['10','20','30','40']\n\n# DEFINE THE MAKEUP FUEL VECTOR FOR THE PROGRESSION\nseed_progression.feed = Eq_Cycle.makeup['seed']\n\n# DEFINE THE STANDARD ORIGEN LIBRARY TO BE USED\nseed_progression.library = 'pwru50'\n\n# LINK PROGRESSION OBJECTS TO KEYS\nEq_Cycle.progressions = {'seed':seed_progression,}\n\n# -----------------------------------------------------------------------------\n# INITIATE EQUILIBRIUM DEPLETION ANALYSIS\n# -----------------------------------------------------------------------------\nEq_Cycle.dBU1 = 1000\nEq_Cycle.dBU2 = 10000\nEq_Cycle.nBU1 = 3\n\nEq_Cycle.search_keff(_BU1_,_BU2_, 1.000)\n"

        #_power_ = '%1.5E' % self.power
        seed_volume = self.dv_seed_entrance + self.dv_seed_expansion + \
            self.dv_seed_active + self.dv_seed_converging + self.dv_seed_defueling
        blanket_volume = self.dv_blanket_entrance + self.dv_blanket_expansion + \
            self.dv_blanket_active + self.dv_blanket_converging + self.dv_blanket_defueling

        ratio = 100*seed_volume/(seed_volume + blanket_volume)

        seed_volume = seed_volume * \
            (self.seed.dv_active/self.seed.V)*(self.seed.TRISO.dv_kernel/self.seed.TRISO.V)

        _power_ = '%1.5E' % (
            self.power_density
            *
            (self.dv_seed_expansion + self.dv_seed_active +
             self.dv_seed_converging) * 1e-6)

        _molar_mass_92235_ = '%1.5E' % (
            seed_volume*self.seed.TRISO.mat_kernel.comp['922350']/4.)
        _molar_mass_92238_ = '%1.5E' % (
            seed_volume*self.seed.TRISO.mat_kernel.comp['922380']/4.)

        _BU1_ = '%1.5E' % (.75*self.seed.TRISO.mat_kernel.comp['922350']/(
            self.seed.TRISO.mat_kernel.comp['922350']+self.seed.TRISO.mat_kernel.comp['922380'])*1.e6)
        _BU2_ = '%1.5E' % (1.25*self.seed.TRISO.mat_kernel.comp['922350']/(
            self.seed.TRISO.mat_kernel.comp['922350']+self.seed.TRISO.mat_kernel.comp['922380'])*1.e6)
        skeleton = skeleton.replace('_power_', _power_)
        skeleton = skeleton.replace('_molar_mass_92235_', _molar_mass_92235_)
        skeleton = skeleton.replace('_molar_mass_92238_', _molar_mass_92238_)
        skeleton = skeleton.replace('_BU1_', _BU1_)
        skeleton = skeleton.replace('_BU2_', _BU2_)

        CHM, density = self.seed.CHM()
        seed_CHM = int(CHM)
        seed_diameter = int(self.seed.TRISO.r_kernel*2.e+04)

        CHM, density = self.blanket.CHM()
        blanket_CHM = int(CHM)

        enrichment = self.seed.TRISO.mat_kernel.comp[
            '922350']/(self.seed.TRISO.mat_kernel.comp['922350']+self.seed.TRISO.mat_kernel.comp['922380'])*100

        dir = '../Equilibrium_Graphite_seed_CHM_%d_dia_%d_enrich_%d_seed_volume_ratio_%2.0f' % (
            seed_CHM, seed_diameter, enrichment, ratio)

        remove = '\\rm %s -r -f'
        # print(remove)
        os.system(remove)

        copy = '\cp -r Equilibrium_BEAU %s' % dir
        # print(copy)
        os.system(copy)

        input_loc = '%s/inp.1' % dir
        open(input_loc, 'w').write(self.input)

        input_loc = '%s/inp' % dir
        open(input_loc, 'w').write(self.input_source)

        mcnp = '#! /bin/sh\n#\n#$ -N Equilibrium_Depletion\n#$ -cwd\n#$ -pe ompi 8\n#$ -S /bin/bash\n#$ -q x.q\n#$ -V\n\nmpiexec -x LD_LIBRARY_PATH  mcnp5.mpi i=inp o=inp.o mc=inp.m runtpe=runtp1\nmv srctp source\nrm inp.o inp.m runtp1\n\n./BEAU7.py\n\n'
        mcnp_loc = '%s/mcnp.sh' % dir
        open(mcnp_loc, 'w').write(mcnp)

        BEAU_loc = '%s/BEAU7.py' % dir
        open(BEAU_loc, 'w').write(skeleton)

        chmod = 'chmod +x %s' % BEAU_loc
        # print(chmod)
        os.system(chmod)


