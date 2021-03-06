#!/usr/bin/env python
from pebble import Pebble
from shield import shield
from heat_exchanger import Heat_Exchanger
import math
import mocup


class PBFHR45:

    def __init__(self):
        self.power = 290  # power of the PBFHR in MWth
        self.mat_coolant = mocup.material()
        self.mat_coolant.comp['30060'] = 20e-6
        self.mat_coolant.comp['30070'] = 2 - 20e-6
        self.mat_coolant.comp['40090'] = 1.
        self.mat_coolant.comp['90190'] = 4.0

        self.coolant_inlet_T = 600. + 273.15
        # coolant inlet temperature (K)
        self.coolant_outlet_T = 700 + 273.15
        # coolant outlet temperature (K)
        self.coolant_temperature_rise = self.coolant_outlet_T - \
            self.coolant_inlet_T
        # temperature rise across the core (C) or (D)
        self.Re = 1200
        # characteristic reynolds number of coolant flow around the pebble

        self.fuel_power_density = 20
        # average power density of active region, converging and expansion
        # regions in MW/m3

        self.wall_PF = 0.6  #.57
        self.bulk_PF = 0.6  #.63

        # Pebble Fuel Object
        self.fuel = Pebble()
        self.fuel.Pebble_PF = self.bulk_PF
        self.fuel.TRISO.dr_kernel = .02
        self.fuel.dr_shell = .1
        self.fuel.generate(300, 1.745)

        # active core dimensions
        self.dr_inner_solid_reflector = 35.
        self.dr_active_fuel = 70.
        self.dr_active_pebble_reflector = 20
        self.dr_active_solid_reflector = 40
        self.dr_core_barrel = 3.
        self.dr_reactor_vessel = 3.
        self.r_reactor_vessel = 175.

        # depeltion parameters
        self.passes = 4
            # number of passes pebbles make through the core
        # radial zoning dimensions
        self.mfp = 3
            # mean free path in (cm)
        self.dt_radial1 = 3.7
            # radial thickness from the inner graphite reflector interms of
            # pebble diameters
        self.df_radial2 = (1./3.)
            # interface between 2nd depletion zone and 3rd depletion zone
            # interms of total radial thickness
        self.dm_radial3 = 3.
            # radial thickness inside the graphite pebble refelctor interms of
            # mean free paths (self.mfp)

        self.pebble_injectors = 3.
        self.dt_buffer = 5.
        # distance from top/bottom of repose piles to the top of the physical
        # pebble region dividers
        self.dz_entrance_shielding = 30.
        # distance from bottom of expansion region to top of the physical pebble
        # region dividers
        self.dt_pebble_reflector = 4.0
        self.dt_defueling_chute = 6.0

        # Heat Exchanger Parameters
        self.DHX_power_density = 1.25
        self.dr_HX_clearance = 0.
        # radial distance from outside of defueling chute to inside of heat
        # exchangers
        self.dz_HX_clearance = 40.
        # axial distance from tope of converging region to the bottom of heat
        # exchangers
        self.decay_heat_fraction = .01
        # ratio of heat source during LOFC or LOHS accident and nominal
        # operation
        self.HX_volume_efficiency = .98
        # cross sectional area fraction of IHX and DHX to total annular cross
        # sectional area
        self.DHX = Heat_Exchanger()
        self.dr_channel = 0.5105
        self.dr_tube = 0.1245
        self.PD = 1.25

        self.converging_angle = math.pi/4.
        self.expansion_angle = math.pi/3.
        self.repose_angle = math.pi/6.
        #
        # ------------
        # \ angle of repose
        #  \  0    /
        #   \     /  < ---- repose pile
        #    \   /
        #     \O/
        #      o
        #    o       < ---- pebbles landing at bottom of pebble bed
        #         o

        self.control_rods = 6
        # number of control rods in inner solid graphite reflector
        self.dt_control_rod_extent = 1.0
        # this is the relative axial extent of the control rods relative to the
        # height of the active region
        self.r_control_rods = 7.0
        # radius of control rod channels (cm)
        self.dt_control_rods = 0.2
        # half thickness of control crucifrom in terms of the control rod
        # channel radius
        self.t_control_rods = 0.8
        # half length of control rod cruciform in terms of the control rod
        # channel radius

        self.control_blades = 12.0
        # number of control blades
        self.dt_control_blades = 0.25
        # half thickness of control blade in terms of a pebble diameter
        self.dm_control_blades = 2.5
        # half length of control blade in long direction in terms of a pebble
        # diameter
        self.dn_control_blades = 0.75
        # half length of control blade in short direction in terms of pebble
        # diameter
        self.dr_control_blades = 4.0
        # clearance for control blades to outer radius of pebble reflector and
        # inner radius of heat exchangers
        self.dr_HX_clearance = self.fuel.r_shell*2*2 * \
            self.dm_control_blades + 2*self.dr_control_blades

        # definition of the shield in the shield object
        self.shield = shield()
        # definition of the b-10 enrichment in w%
        self.shield.setEnrichment(.5)
        # radial thickness of the shield
        self.shield_dr = 3.

        # definition of the porous reflector regions
        #	This is the region composed of coolant and graphite reflector
        # radial thickness of the porous reflector region (cm)
        self.inner_porous_reflector_dr = 10.
        self.outer_porous_reflector_dr = 10.
        # volume fraction of the porous reflector
        #	1 - volume fraction of coolant in
        self.inner_porous_reflector_PF = .6
        self.outer_porous_reflector_PF = .6

        # definition of bounds of the inner porous reflector region nad outer porous reflector region
        # zf_inner_porous is the fraction of the height (from the bottom) of the active region that the inner porous region occupies
        #
        #		|       |         |      |
        #       |       |		  | 8021_|
        #       |_      | Pebbles |    |
        #		  | 8020|         |    |
        # inner   |     |         |    | outer
        # refl.   |     |         |    | refl.
        #
        # 	ie 0 corrisponds to a short inlet and 1 corrisponds to a tall inlet
        self.zf_inner_porous = 1 - 1e-6
        # zf_inner_porous is the fraction of the height (from the top) of the
        # active region that the outer porous region occupies
        self.zf_outer_porous = 1 - 1e-6
        # ie 0 corrisponds to a shallow outlet and 1 corrisponds to a deep
        # outlet

        self.tallies = []
        # enter
        # Dose or dose for dose mesh tallies in the inner graphite reflector, otuer graphite reflector, heat exchanger, core barrel and reactor vessel
        # Spectrum or spectrum for a neutron spectrum tally in the fuel
        # Mocup or mocup for the mocup tallies and single isotope materials
        # Leakage or leakage for the tallies to calculate leakage from the fuel pebble region
        # Economy or economy for the tallies to calculate the neutron economics
        # and FHR 4 factor formula

    def generate(self):
        self.fuel.volume()
        self.dr_HX_clearance = self.fuel.r_shell*2*2 * \
            self.dm_control_blades + 2*self.dr_control_blades

        core_volume = self.power / self.fuel_power_density * 1e+6
        print 'core_volume%f' %core_volume
        self.r_inner_solid_reflector = self.dr_inner_solid_reflector
        self.r_active_fuel = self.r_inner_solid_reflector + self.dr_active_fuel
        self.r_active_pebble_reflector = self.r_active_fuel + \
            self.dr_active_pebble_reflector
        fuel_cross_section_fraction = (
            self.r_active_fuel**2 - self.r_inner_solid_reflector**2.)/(
            self.r_active_pebble_reflector**2. - self.r_inner_solid_reflector**2.)
        self.r_active_solid_reflector = self.r_active_pebble_reflector + \
            self.dr_active_solid_reflector
        self.r_core_barrel = self.r_active_solid_reflector + self.dr_core_barrel
        self.r_active_fuel1 = self.r_inner_solid_reflector + \
            self.dt_radial1*self.fuel.r_shell*2.
        self.r_active_fuel2 = self.r_inner_solid_reflector + \
            self.dr_active_fuel*self.df_radial2
        self.r_active_fuel3 = self.r_inner_solid_reflector + \
            self.dr_active_fuel - self.dm_radial3*self.mfp
        self.r_active_fuel4 = self.r_inner_solid_reflector + self.dr_active_fuel
        self.r_active_pebble_reflector1 = self.r_active_fuel4 + \
            self.dr_active_pebble_reflector - self.dt_radial1*self.fuel.r_shell*2.
        self.r_active_pebble_reflector2 = self.r_active_fuel4 + \
            self.dr_active_pebble_reflector

        # determine effective area fractions for entrance and defueling chutes - where there is uniform porosity
        # active height, volume of pebbles are constant so these are terms are
        # dropped
        self.a1 = math.pi * \
            (self.r_active_fuel1**2. - self.r_inner_solid_reflector**2.)*self.wall_PF
        self.a2 = math.pi * \
            (self.r_active_fuel2**2. - self.r_active_fuel1**2.)*self.bulk_PF
        self.a3 = math.pi * \
            (self.r_active_fuel3**2. - self.r_active_fuel2**2.)*self.bulk_PF
        self.a4 = math.pi * \
            (self.r_active_fuel4**2. - self.r_active_fuel3**2.)*self.bulk_PF
        self.a5 = math.pi * \
            (self.r_active_pebble_reflector1**2. - self.r_active_fuel4**2.)*self.bulk_PF
        self.a6 = math.pi*(self.r_active_pebble_reflector2 **
                           2. - self.r_active_pebble_reflector1**2.)*self.wall_PF

        self.df_a1 = self.a1 / \
            (self.a1 + self.a2 + self.a3 + self.a4 + self.a5 + self.a6)
        self.df_a2 = self.a2 / \
            (self.a1 + self.a2 + self.a3 + self.a4 + self.a5 + self.a6)
        self.df_a3 = self.a3 / \
            (self.a1 + self.a2 + self.a3 + self.a4 + self.a5 + self.a6)
        self.df_a4 = self.a4 / \
            (self.a1 + self.a2 + self.a3 + self.a4 + self.a5 + self.a6)
        self.df_a5 = (self.a5 + self.a6) / \
            (self.a1 + self.a2 + self.a3 + self.a4 + self.a5 + self.a6)

        fuel_cross_section_fraction = self.df_a1 + \
            self.df_a2 + self.df_a3 + self.df_a4

        # -------------------------------------------------------------------------------------
        # solve for pebble reflector radius in entrance region using quadratic formula
        # quadradic coefficients
        # -------------------------------------------------------------------------------------

        # 0 = a*r**2 + b*r + c
        a = 1. - fuel_cross_section_fraction
        b = -2.*(self.fuel.r_shell*2)*self.dt_pebble_reflector
        c = (self.fuel.r_shell*2)**2*self.dt_pebble_reflector**2 + \
            (fuel_cross_section_fraction-1)*self.r_inner_solid_reflector**2.

        self.r_entrance_pebble_reflector = (-b + (b**2. - 4.*a*c)**0.5)/(2.*a)

        self.r_entrance_fuel = self.r_entrance_pebble_reflector - \
            self.dt_pebble_reflector*self.fuel.r_shell*2
        self.dz_expansion = 1./math.tan(math.pi/2. - self.expansion_angle) * \
            (self.r_active_pebble_reflector - self.r_entrance_pebble_reflector)
        self.r_defueling_pebble_reflector = self.r_inner_solid_reflector + \
            self.dt_defueling_chute*(self.fuel.r_shell*2)
        self.r_defueling_fuel = (
            (self.r_defueling_pebble_reflector**2 -
             self.r_inner_solid_reflector**2) *
            fuel_cross_section_fraction +
            self.r_inner_solid_reflector**2)**0.5
        self.dz_converging = 1./math.tan(math.pi/2. - self.converging_angle) * \
            (self.r_active_pebble_reflector - self.r_defueling_pebble_reflector)

        self.dv_expansion = (1./3.)*math.pi*(self.dz_expansion/(self.r_active_fuel-self.r_entrance_fuel)) * \
            (self.r_active_fuel**3. - self.r_entrance_fuel**3.) - math.pi*self.r_inner_solid_reflector**2.*self.dz_expansion
        self.dv_converging = (1./3.)*math.pi*(self.dz_converging/(self.r_active_fuel-self.r_defueling_fuel)) * \
            (self.r_active_fuel**3. - self.r_defueling_fuel**3.) - math.pi*self.r_inner_solid_reflector**2.*self.dz_converging

        self.dv_active = core_volume - self.dv_expansion - self.dv_converging
        self.dz_active = self.dv_active / \
            (math.pi*(self.r_active_fuel**2. - self.r_inner_solid_reflector**2.))

        # -------------------------------------------------------------------------------------
        # entrance
        # -------------------------------------------------------------------------------------

        self.dz_pile = self.r_entrance_pebble_reflector*math.pi / \
            self.pebble_injectors*math.tan(self.repose_angle)
        dz_buffer = 2*self.fuel.r_shell*self.dt_buffer
        self.dz_entrance = self.dz_pile + dz_buffer + self.dz_entrance_shielding

        area = math.pi * \
            (self.r_entrance_pebble_reflector**2 - self.r_inner_solid_reflector**2)
        a1 = math.pi*self.r_inner_solid_reflector**2 + self.df_a1*area
        a2 = math.pi*self.r_inner_solid_reflector**2 + \
            (self.df_a2 + self.df_a1)*area
        a3 = math.pi*self.r_inner_solid_reflector**2 + \
            (self.df_a3 + self.df_a2 + self.df_a1)*area
        a4 = math.pi*self.r_inner_solid_reflector**2 + \
            (self.df_a4 + self.df_a3 + self.df_a2 + self.df_a1)*area

        self.r_entrance_fuel1 = ((a1)/math.pi)**0.5
        self.r_entrance_fuel2 = ((a2)/math.pi)**0.5
        self.r_entrance_fuel3 = ((a3)/math.pi)**0.5
        self.r_entrance_fuel4 = ((a4)/math.pi)**0.5

        # -------------------------------------------------------------------------------------
        # defueling chute
        # -------------------------------------------------------------------------------------

        area = math.pi * \
            (self.r_defueling_pebble_reflector**2. - self.r_inner_solid_reflector**2)
        a1 = math.pi*self.r_inner_solid_reflector**2 + self.df_a1*area
        a2 = math.pi*self.r_inner_solid_reflector**2 + \
            (self.df_a2 + self.df_a1)*area
        a3 = math.pi*self.r_inner_solid_reflector**2 + \
            (self.df_a3 + self.df_a2 + self.df_a1)*area
        a4 = math.pi*self.r_inner_solid_reflector**2 + \
            (self.df_a4 + self.df_a3 + self.df_a2 + self.df_a1)*area

        self.r_defueling_fuel1 = ((a1)/math.pi)**0.5
        self.r_defueling_fuel2 = ((a2)/math.pi)**0.5
        self.r_defueling_fuel3 = ((a3)/math.pi)**0.5
        self.r_defueling_fuel4 = ((a4)/math.pi)**0.5

        # -------------------------------------------------------------------------------------
        # heat exchangers
        # -------------------------------------------------------------------------------------

        HX_volume = (
            self.power *
            1e6 *
            self.decay_heat_fraction /
            self.DHX_power_density)
        self.r_defueling_reflector = self.r_defueling_pebble_reflector + \
            self.dr_HX_clearance
        self.r_HX = self.r_reactor_vessel - self.dr_reactor_vessel
        if self.r_core_barrel > self.r_HX:
            print('!!!Geometry Error!!!\n!!!Core Barrel Reactor Vessel Overlap!!!')
            poop
        HX_cross_section = math.pi * \
            (self.r_HX**2. - self.r_defueling_reflector**2.)

        partial_core_volume = core_volume + math.pi * \
            (self.r_entrance_pebble_reflector**2 - self.r_inner_solid_reflector**2)*self.dz_entrance

        self.dz_HX = (
            partial_core_volume + HX_volume + math.pi
            *
            (self.r_defueling_pebble_reflector ** 2. -
             self.r_inner_solid_reflector ** 2) * self.dz_HX_clearance) / (
            (HX_cross_section * self.HX_volume_efficiency - math.pi
             *
             (self.r_defueling_pebble_reflector ** 2. -
              self.r_inner_solid_reflector ** 2)))

        self.dv_DHX = HX_volume/(HX_cross_section*self.dz_HX)
        self.dv_decay_storage = (
            partial_core_volume + math.pi
            *
            (self.r_defueling_pebble_reflector ** 2. -
             self.r_inner_solid_reflector ** 2)
            * (self.dz_HX_clearance + self.dz_HX)) / (
            HX_cross_section * self.dz_HX)
        self.dv_reflector = (1. - self.dv_DHX - self.dv_decay_storage)

        self.mat_reflector = mocup.material()
        self.mat_reflector.comp['60120'] = 1.74/12.

        self.mat_HX = self.DHX.mat()*self.dv_DHX + self.fuel.mat_pebble * \
            self.dv_decay_storage + self.mat_reflector*self.dv_reflector

        self.Z = self.dz_entrance + self.dz_expansion + self.dz_active + \
            self.dz_converging + self.dz_HX_clearance + self.dz_HX
        # print(self.Z)

    def MCNP5(self):
        coolant_temperature = (self.coolant_inlet_T + self.coolant_outlet_T)/2.
        #self.fuel.Pebble_PF = self.bulk_PF
        self.fuel.HT(coolant_temperature, self.Re, self.fuel_power_density)
        self.fuel.TRISO.matrix()
        _fuel_density_ = '%1.5E' % self.fuel.TRISO.mat_kernel.mass()
        _fuel_temperature_ = '%1.5E' % (self.fuel.T_kernel*8.6173324e-11)

        # calculate fuel volumes
        # entrance
        _fuel_volume11_ = '%1.5E' % (
            math.pi * (self.dz_entrance - self.dz_pile / 2.)
            * (self.r_entrance_fuel1 ** 2 - self.r_inner_solid_reflector ** 2) *
            self.wall_PF * (self.fuel.r_active ** 3. - self.fuel.r_core ** 3.) /
            self.fuel.r_shell ** 3.
            * (self.fuel.TRISO.dv_kernel / self.fuel.TRISO.V) / self.passes)
        _fuel_volume21_ = '%1.5E' % (
            math.pi * (self.dz_entrance - self.dz_pile / 2.)
            * (self.r_entrance_fuel2 ** 2 - self.r_entrance_fuel1 ** 2) *
            self.wall_PF * (self.fuel.r_active ** 3. - self.fuel.r_core ** 3.) /
            self.fuel.r_shell ** 3.
            * (self.fuel.TRISO.dv_kernel / self.fuel.TRISO.V) / self.passes)
        _fuel_volume31_ = '%1.5E' % (
            math.pi * (self.dz_entrance - self.dz_pile / 2.)
            * (self.r_entrance_fuel3 ** 2 - self.r_entrance_fuel2 ** 2) *
            self.wall_PF * (self.fuel.r_active ** 3. - self.fuel.r_core ** 3.) /
            self.fuel.r_shell ** 3.
            * (self.fuel.TRISO.dv_kernel / self.fuel.TRISO.V) / self.passes)
        _fuel_volume41_ = '%1.5E' % (
            math.pi * (self.dz_entrance - self.dz_pile / 2.)
            * (self.r_entrance_fuel4 ** 2 - self.r_entrance_fuel3 ** 2) *
            self.wall_PF * (self.fuel.r_active ** 3. - self.fuel.r_core ** 3.) /
            self.fuel.r_shell ** 3.
            * (self.fuel.TRISO.dv_kernel / self.fuel.TRISO.V) / self.passes)
        # expansion
        vol0 = math.pi*self.dz_expansion*self.r_inner_solid_reflector**2
        vol1 = (1. / 3.) * math.pi * self.dz_expansion / (self.r_active_fuel1 -
                                                          self.r_entrance_fuel1) * (self.r_active_fuel1 ** 3. -
                                                                                    self.r_entrance_fuel1 ** 3.)
        vol2 = (1. / 3.) * math.pi * self.dz_expansion / (self.r_active_fuel2 -
                                                          self.r_entrance_fuel2) * (self.r_active_fuel2 ** 3. -
                                                                                    self.r_entrance_fuel2 ** 3.)
        vol3 = (1. / 3.) * math.pi * self.dz_expansion / (self.r_active_fuel3 -
                                                          self.r_entrance_fuel3) * (self.r_active_fuel3 ** 3. -
                                                                                    self.r_entrance_fuel3 ** 3.)
        vol4 = (1. / 3.) * math.pi * self.dz_expansion / (self.r_active_fuel4 -
                                                          self.r_entrance_fuel4) * (self.r_active_fuel4 ** 3. -
                                                                                    self.r_entrance_fuel4 ** 3.)
        _fuel_volume12_ = '%1.5E' % (
            (vol1 - vol0) * self.wall_PF
            * (self.fuel.r_active ** 3. - self.fuel.r_core ** 3.) /
            self.fuel.r_shell ** 3.
            * (self.fuel.TRISO.dv_kernel / self.fuel.TRISO.V) / self.passes)
        _fuel_volume22_ = '%1.5E' % (
            (vol2 - vol1) * self.wall_PF
            * (self.fuel.r_active ** 3. - self.fuel.r_core ** 3.) /
            self.fuel.r_shell ** 3.
            * (self.fuel.TRISO.dv_kernel / self.fuel.TRISO.V) / self.passes)
        _fuel_volume32_ = '%1.5E' % (
            (vol3 - vol2) * self.wall_PF
            * (self.fuel.r_active ** 3. - self.fuel.r_core ** 3.) /
            self.fuel.r_shell ** 3.
            * (self.fuel.TRISO.dv_kernel / self.fuel.TRISO.V) / self.passes)
        _fuel_volume42_ = '%1.5E' % (
            (vol4 - vol3) * self.wall_PF
            * (self.fuel.r_active ** 3. - self.fuel.r_core ** 3.) /
            self.fuel.r_shell ** 3.
            * (self.fuel.TRISO.dv_kernel / self.fuel.TRISO.V) / self.passes)

        # active

        _fuel_volume13_ = '%1.5E' % (
            math.pi * self.dz_active
            * (self.r_active_fuel1 ** 2 - self.r_inner_solid_reflector ** 2) *
            self.wall_PF * (self.fuel.r_active ** 3. - self.fuel.r_core ** 3.) /
            self.fuel.r_shell ** 3.
            * (self.fuel.TRISO.dv_kernel / self.fuel.TRISO.V) / self.passes)
        _fuel_volume23_ = '%1.5E' % (
            math.pi * self.dz_active
            * (self.r_active_fuel2 ** 2 - self.r_active_fuel1 ** 2) *
            self.bulk_PF * (self.fuel.r_active ** 3. - self.fuel.r_core ** 3.) /
            self.fuel.r_shell ** 3.
            * (self.fuel.TRISO.dv_kernel / self.fuel.TRISO.V) / self.passes)
        _fuel_volume33_ = '%1.5E' % (
            math.pi * self.dz_active
            * (self.r_active_fuel3 ** 2 - self.r_active_fuel2 ** 2) *
            self.bulk_PF * (self.fuel.r_active ** 3. - self.fuel.r_core ** 3.) /
            self.fuel.r_shell ** 3.
            * (self.fuel.TRISO.dv_kernel / self.fuel.TRISO.V) / self.passes)
        _fuel_volume43_ = '%1.5E' % (
            math.pi * self.dz_active
            * (self.r_active_fuel4 ** 2 - self.r_active_fuel3 ** 2) *
            self.bulk_PF * (self.fuel.r_active ** 3. - self.fuel.r_core ** 3.) /
            self.fuel.r_shell ** 3.
            * (self.fuel.TRISO.dv_kernel / self.fuel.TRISO.V) / self.passes)

        # converging

        vol0 = math.pi*self.dz_converging*self.r_inner_solid_reflector**2
        vol1 = (1. / 3.) * math.pi * self.dz_converging / \
            (self.r_active_fuel1 - self.r_defueling_fuel1) *\
            (self.r_active_fuel1 ** 3. - self.r_defueling_fuel1 ** 3.)
        vol2 = (1. / 3.) * math.pi * self.dz_converging /\
            (self.r_active_fuel2 - self.r_defueling_fuel2) *\
            (self.r_active_fuel2 ** 3. - self.r_defueling_fuel2 ** 3.)
        vol3 = (1. / 3.) * math.pi * self.dz_converging / \
            (self.r_active_fuel3 -  self.r_defueling_fuel3) * \
            (self.r_active_fuel3 ** 3. - self.r_defueling_fuel3 ** 3.)
        vol4 = (1. / 3.) * math.pi * self.dz_converging /\
            (self.r_active_fuel4 - self.r_defueling_fuel4) *\
            (self.r_active_fuel4 ** 3. - self.r_defueling_fuel4 ** 3.)

        _fuel_volume14_ = '%1.5E' % (
            (vol1 - vol0) * self.wall_PF
            * (self.fuel.r_active ** 3. - self.fuel.r_core ** 3.) /
            self.fuel.r_shell ** 3.
            * (self.fuel.TRISO.dv_kernel / self.fuel.TRISO.V) / self.passes)
        _fuel_volume24_ = '%1.5E' % (
            (vol2 - vol1) * self.wall_PF
            * (self.fuel.r_active ** 3. - self.fuel.r_core ** 3.) /
            self.fuel.r_shell ** 3.
            * (self.fuel.TRISO.dv_kernel / self.fuel.TRISO.V) / self.passes)
        _fuel_volume34_ = '%1.5E' % (
            (vol3 - vol2) * self.wall_PF
            * (self.fuel.r_active ** 3. - self.fuel.r_core ** 3.) /
            self.fuel.r_shell ** 3.
            * (self.fuel.TRISO.dv_kernel / self.fuel.TRISO.V) / self.passes)
        _fuel_volume44_ = '%1.5E' % (
            (vol4 - vol3) * self.wall_PF
            * (self.fuel.r_active ** 3. - self.fuel.r_core ** 3.) /
            self.fuel.r_shell ** 3.
            * (self.fuel.TRISO.dv_kernel / self.fuel.TRISO.V) / self.passes)

        # defueling

        _fuel_volume15_ = '%1.5E' % (
            math.pi * (self.dz_HX_clearance + self.dz_HX)
            *
            (self.r_defueling_fuel1 ** 2 - self.r_inner_solid_reflector ** 2) *
            self.wall_PF * (self.fuel.r_active ** 3. - self.fuel.r_core ** 3.) /
            self.fuel.r_shell ** 3.
            * (self.fuel.TRISO.dv_kernel / self.fuel.TRISO.V) / self.passes)
        _fuel_volume25_ = '%1.5E' % (
            math.pi * (self.dz_HX_clearance + self.dz_HX)
            * (self.r_defueling_fuel2 ** 2 - self.r_defueling_fuel1 ** 2) *
            self.wall_PF * (self.fuel.r_active ** 3. - self.fuel.r_core ** 3.) /
            self.fuel.r_shell ** 3.
            * (self.fuel.TRISO.dv_kernel / self.fuel.TRISO.V) / self.passes)
        _fuel_volume35_ = '%1.5E' % (
            math.pi * (self.dz_HX_clearance + self.dz_HX)
            * (self.r_defueling_fuel3 ** 2 - self.r_defueling_fuel2 ** 2) *
            self.wall_PF * (self.fuel.r_active ** 3. - self.fuel.r_core ** 3.) /
            self.fuel.r_shell ** 3.
            * (self.fuel.TRISO.dv_kernel / self.fuel.TRISO.V) / self.passes)
        _fuel_volume45_ = '%1.5E' % (
            math.pi * (self.dz_HX_clearance + self.dz_HX)
            * (self.r_defueling_fuel4 ** 2 - self.r_defueling_fuel3 ** 2) *
            self.wall_PF * (self.fuel.r_active ** 3. - self.fuel.r_core ** 3.) /
            self.fuel.r_shell ** 3.
            * (self.fuel.TRISO.dv_kernel / self.fuel.TRISO.V) / self.passes)

        _matrix_density_ = '%1.5E' % self.fuel.TRISO.mat_hmatrix.mass()
        _matrix_temperature_ = '%1.5E' % (self.fuel.T_active*8.6173324e-11)

        _core_density_ = '%1.5E' % self.fuel.mat_core.mass()
        _shell_density_ = '%1.5E' % self.fuel.mat_shell.mass()
        _core_temperature_ = '%1.5E' % (self.fuel.T_core*8.6173324e-11)
        _shell_temperature_ = '%1.5E' % (self.fuel.T_shell*8.6173324e-11)

        _coolant_outlet_density_ = '%1.5E' % (
            2.28 - 0.000488*(self.coolant_outlet_T - 273.15))
        _coolant_density_ = '%1.5E' % (
            2.28 - 0.000488*(coolant_temperature - 273.15))
        _coolant_inlet_density_ = '%1.5E' % (
            2.28 - 0.000488*(self.coolant_inlet_T - 273.15))

        _coolant_outlet_temperature_ = '%1.5E' % (
            self.coolant_outlet_T*8.6173324e-11)
        _coolant_temperature_ = '%1.5E' % ((coolant_temperature*8.6173324e-11))
        _coolant_inlet_temperature_ = '%1.5E' % (
            self.coolant_inlet_T*8.6173324e-11)

        _kernel_radius_ = '%1.5E' % self.fuel.TRISO.r_kernel
        _triso_hpitch_ = '%1.5E' % self.fuel.TRISO.hpitch

        self.fuel.Pebble_PF = self.wall_PF
        self.fuel.volume()
        _pebble_hpitch1_ = '%1.5E' % self.fuel.h_pitch
        self.fuel.Pebble_PF = self.bulk_PF
        self.fuel.volume()
        _pebble_hpitch2_ = '%1.5E' % self.fuel.h_pitch
        _pebble_core_radius_ = '%1.5E' % self.fuel.r_core
        _pebble_active_radius_ = '%1.5E' % self.fuel.r_active
        _pebble_shell_radius_ = '%1.5E' % self.fuel.r_shell

        _inner_solid_reflector_radius_ = '%1.5E' % self.r_inner_solid_reflector
        _inner_solid_reflector_interface_ = '%1.5E' % (
            self.
            r_inner_solid_reflector -
            self.
            inner_porous_reflector_dr)
        _outer_solid_refelctor_radius_ = '%1.5E' % self.r_active_solid_reflector
        _outer_solid_refelctor_interface_ = '%1.5E' % (
            self.r_active_pebble_reflector + self.outer_porous_reflector_dr)
        _shield_radius_ = '%1.5E' % (
            self.r_active_solid_reflector - self.shield_dr)
        _core_barrel_radius_ = '%1.5E' % self.r_core_barrel
        _down_comer_radius_ = '%1.5E' % self.r_HX
        _reactor_vessel_radius_ = '%1.5E' % self.r_reactor_vessel

        _entrance_fuel_radius1_ = '%1.5E' % self.r_entrance_fuel1
        _entrance_fuel_radius2_ = '%1.5E' % self.r_entrance_fuel2
        _entrance_fuel_radius3_ = '%1.5E' % self.r_entrance_fuel3
        _entrance_fuel_radius4_ = '%1.5E' % self.r_entrance_fuel4
        _entrance_pebble_reflector_radius_ = '%1.5E' % self.r_entrance_pebble_reflector

        t = (self.r_active_fuel1 - self.r_entrance_fuel1)/self.dz_expansion
        z = (self.dz_entrance + self.dz_expansion) - self.r_active_fuel1/t
        _expansion_fuel_cone1_ = '%1.5E %1.5E 1' % (z, t**2)
        t = (self.r_active_fuel2 - self.r_entrance_fuel2)/self.dz_expansion
        z = (self.dz_entrance + self.dz_expansion) - self.r_active_fuel2/t
        _expansion_fuel_cone2_ = '%1.5E %1.5E 1' % (z, t**2)
        t = (self.r_active_fuel3 - self.r_entrance_fuel3)/self.dz_expansion
        z = (self.dz_entrance + self.dz_expansion) - self.r_active_fuel3/t
        _expansion_fuel_cone3_ = '%1.5E %1.5E 1' % (z, t**2)
        t = (self.r_active_fuel4 - self.r_entrance_fuel4)/self.dz_expansion
        z = (self.dz_entrance + self.dz_expansion) - self.r_active_fuel4/t
        _expansion_fuel_cone4_ = '%1.5E %1.5E 1' % (z, t**2)
        t = (self.r_active_pebble_reflector2 -
             self.r_entrance_pebble_reflector)/self.dz_expansion
        z = (self.dz_entrance + self.dz_expansion) - \
            self.r_active_pebble_reflector2/t
        _expansion_pebble_reflector_core_ = '%1.5E %1.5E 1' % (z, t**2)

        _active_fuel_radius1_ = '%1.5E' % self.r_active_fuel1
        _active_fuel_radius2_ = '%1.5E' % self.r_active_fuel2
        _active_fuel_radius3_ = '%1.5E' % self.r_active_fuel3
        _active_fuel_radius4_ = '%1.5E' % self.r_active_fuel4
        _active_pebble_reflector_radius1_ = '%1.5E' % self.r_active_pebble_reflector1
        _active_pebble_reflector_radius2_ = '%1.5E' % self.r_active_pebble_reflector2

        t = (self.r_active_fuel1 - self.r_defueling_fuel1)/self.dz_converging
        z = (self.dz_active + self.dz_entrance + self.dz_expansion) + \
            self.r_active_fuel1/t
        _converging_fuel_cone1_ = '%1.5E %1.5E -1' % (z, t**2)
        t = (self.r_active_fuel2 - self.r_defueling_fuel2)/self.dz_converging
        z = (self.dz_active + self.dz_entrance + self.dz_expansion) + \
            self.r_active_fuel2/t
        _converging_fuel_cone2_ = '%1.5E %1.5E -1' % (z, t**2)
        t = (self.r_active_fuel3 - self.r_defueling_fuel3)/self.dz_converging
        z = (self.dz_active + self.dz_entrance + self.dz_expansion) + \
            self.r_active_fuel3/t
        _converging_fuel_cone3_ = '%1.5E %1.5E -1' % (z, t**2)
        t = (self.r_active_fuel4 - self.r_defueling_fuel4)/self.dz_converging
        z = (self.dz_active + self.dz_entrance + self.dz_expansion) + \
            self.r_active_fuel4/t
        _converging_fuel_cone4_ = '%1.5E %1.5E -1' % (z, t**2)
        t = (self.r_active_pebble_reflector -
             self.r_defueling_pebble_reflector)/self.dz_converging
        z = (self.dz_active + self.dz_entrance + self.dz_expansion) + \
            self.r_active_pebble_reflector/t
        _converging_pebble_reflector_cone_ = '%1.5E %1.5E -1' % (z, t**2)

        _defueling_fuel_radius1_ = '%1.5E' % self.r_defueling_fuel1
        _defueling_fuel_radius2_ = '%1.5E' % self.r_defueling_fuel2
        _defueling_fuel_radius3_ = '%1.5E' % self.r_defueling_fuel3
        _defueling_fuel_radius4_ = '%1.5E' % self.r_defueling_fuel4
        _defueling_pebble_reflector_radius_ = '%1.5E' % self.r_defueling_pebble_reflector

        if (self.dz_entrance/self.r_core_barrel) > 0.5:
            #
            # -------------------------
            # |                        |
            # |                        |
            #  \                      /
            #    \                  /
            #      \              /
            #         \        /
            #              -

            z = self.r_core_barrel*0.5

        else:

            #
            # -------------------------
            #  \                      /
            #    \                  /
            #      \              /
            #         \        /
            #              -

            z = self.dz_entrance

        _entrance2_axial_height_ = '%1.5E' % z

        A = 1./(self.r_core_barrel**2.)
        B = 1./(self.r_core_barrel**2.)
        C = (1. - A*self.r_entrance_pebble_reflector**2.)/z**2.

        #print('A: %1.5E' % A)
        #print('C: %1.5E' % C)
        #print('r entrance pebble reflector %1.5E' % self.r_entrance_pebble_reflector)
        #print('r core barrel %1.5E' % self.r_core_barrel)
        # print(z)
        # print(self.dz_entrance)

        Z = (1./C)**0.5

        _entrance_bottom_ = '%1.5E %1.5E %1.5E\n           0 0 0 -1 0 0 %1.5E' % (
            A, B, C, z)

        Z += (self.r_HX - self.r_core_barrel)

        A = 1./(self.r_HX**2.)
        B = 1./(self.r_HX**2.)
        C = 1./(Z**2.)

        _plenum_coolant_surface_ = '''
        %1.5E %1.5E %1.5E\n
        0 0 0 -1 0 0 %1.5E' % (
            A, B, C, z)
        '''
        Z += self.dr_reactor_vessel
        A = 1./(self.r_reactor_vessel**2.)
        B = 1./(self.r_reactor_vessel**2.)
        C = 1/(Z**2.)

        _plenum_reactor_vessel_ = '%1.5E %1.5E %1.5E\n           0 0 0 -1 0 0 %1.5E' % (
            A, B, C, z)

        self.Z += (Z - z)

        _repose_height_ = '%1.5E' % (self.dz_pile/2.)
        _entrance_axial_height_ = '%1.5E' % (self.dz_entrance)
        _expansion_axial_height_ = '%1.5E' % (
            self.dz_entrance + self.dz_expansion)
        _active_axial_height_ = '%1.5E' % (
            self.dz_entrance + self.dz_expansion + self.dz_active)
        _bottom_control_rod_ = '%1.5E' % (
            self.dz_entrance + self.dz_expansion + self.dz_active
            * (1. - self.dt_control_rod_extent) / 2.)
        _top_control_rod_ = '%1.5E' % (
            self.dz_entrance + self.dz_expansion + self.dz_active
            * (1. + self.dt_control_rod_extent) / 2.)
        _top_porous_inner_reflector_ = '%1.5E' % (
            self.dz_entrance + (self.dz_expansion + self.dz_active + self.dz_converging)*self.zf_inner_porous)
        _bottom_porous_outer_reflector_ = '%1.5E' % (
            self.dz_entrance + self.dz_expansion + self.dz_active
            * (1. - self.zf_inner_porous))
        _converging_axial_height_ = '%1.5E' % (
            self.dz_entrance + self.dz_expansion + self.dz_active +
            self.dz_converging)
        _defueling_axial_height_ = '%1.5E' % (self.Z)
        self.total_height = self.Z + (Z - z)


        mat = mocup.material()
        # mat.addnux()
        fuel = (self.fuel.TRISO.mat_kernel*.602214129 + mat)

        _fuel_materials_ = ''
        _fuel_cells_ = ''
        _matrix_cells_ = ''
        _pebble_shell_cells_ = ''
        _pebble_core_cells_ = ''
        _coolant_cells_ = ''
        _fuel_tally_ = ''

        for z in range(1, 6):
            for r in range(1, 5):
                _RZ_ = '%d%d' % (r, z)
                _fuel_cells_ += '          '
                _matrix_cells_ += '          '
                for b in range(1, self.passes + 1):
                    m = int(r*1e4+z*1e3+b*1e2)
                    _fuel_materials_ += fuel.mcf(m,
                                                 tmp=self.fuel.T_kernel) + '\n'
                    _fuel_cells_ += '%d ' % m
                    _matrix_cells_ += '%d ' % (m+1)
                    _fuel_tally_ += '          (-1 %d (-7 -6) (-6:-2))\n' % m

                _fuel_cells_ += '\n'
                _matrix_cells_ += '\n'
                _pebble_core_cells_ += ' '*10 + '_RZ_00 _RZ_03 _RZ_06 _RZ_09\n'
                _pebble_core_cells_ += ' '*10 + '_RZ_20 _RZ_23 _RZ_26 _RZ_29\n'
                _pebble_core_cells_ += ' '*10 + \
                    '_RZ_40 _RZ_03 _RZ_50 _RZ_53 _RZ_60 _RZ_63\n'
                _pebble_core_cells_ = _pebble_core_cells_.replace('_RZ_', _RZ_)

                _pebble_shell_cells_ += ' '*10 + '_RZ_02 _RZ_05 _RZ_08 _RZ_11\n'
                _pebble_shell_cells_ += ' '*10 + '_RZ_22 _RZ_25 _RZ_28 _RZ_31\n'
                _pebble_shell_cells_ += ' '*10 + \
                    '_RZ_42 _RZ_45 _RZ_52 _RZ_55 _RZ_62 _RZ_65\n'
                _pebble_shell_cells_ = _pebble_shell_cells_.replace(
                    '_RZ_',
                    _RZ_)

                _coolant_cells_ += ' '*10 + _RZ_ + '80\n'

        _heat_exchanger_radius_ = '%1.5E' % (
            self.r_defueling_pebble_reflector + self.fuel.r_shell * 2 * 2 *
            self.dm_control_blades + 2 * self.dr_control_blades)
        _heat_exchanger_axial_height_ = '%1.5E' % (
            self.dz_entrance + self.dz_expansion + self.dz_active +
            self.dz_converging + self.dz_HX_clearance)
        _heat_exchanger_density_ = '%1.5E' % self.mat_HX.mass()

        graphite = mocup.material()
        graphite.comp['60120'] = 1.74/12.

        flibe = self.mat_coolant

        control = mocup.material()
        control.comp['50100'] = .199
        control.comp['50110'] = .801
        control.comp['60120'] = 4.
        control = control*(2.4/control.mass())

        inlet_porous_reflector = graphite*self.inner_porous_reflector_PF + flibe * \
            (1 - self.inner_porous_reflector_PF)*((2.28 - 0.000488*(self.coolant_inlet_T-273.15))/flibe.mass())
        outlet_porous_reflector = graphite*self.inner_porous_reflector_PF + flibe * \
            (1 - self.inner_porous_reflector_PF)*((2.28 - 0.000488*(self.coolant_outlet_T-273.15))/flibe.mass())

        _density_inner_porous_reflector_ = '%1.5E' % inlet_porous_reflector.mass(
            )
        _density_outer_porous_reflector_ = '%1.5E' % outlet_porous_reflector.mass(
            )
        _inner_solid_reflector_density_ = '%1.5E' % graphite.mass()

        _graphite_inlet_ = graphite.mcf(
            1,
            tmp=self.coolant_inlet_T,
            scat='grph')
        _graphite_outlet_ = graphite.mcf(
            2,
            tmp=self.coolant_outlet_T,
            scat='grph')
        _porous_graphite_inlet_ = inlet_porous_reflector.mcf(
            5,
            tmp=self.coolant_inlet_T,
            scat='grph')
        _porous_graphite_outlet_ = outlet_porous_reflector.mcf(
            6,
            tmp=self.coolant_outlet_T,
            scat='grph')
        _graphite_pebble_core_ = graphite.mcf(
            3,
            tmp=self.fuel.T_core,
            scat='grph')
        _graphite_pebble_shell_ = graphite.mcf(
            4,
            tmp=self.fuel.T_shell,
            scat='grph')
        _triso_matrix_ = self.fuel.TRISO.mat_hmatrix.mcf(
            11,
            tmp=self.fuel.T_active,
            scat='grph')
        _inlet_coolant_material_ = flibe.mcf(23, tmp=self.coolant_inlet_T)
        _coolant_material_ = flibe.mcf(
            24,
            tmp=(
                self.coolant_inlet_T +
                self.coolant_outlet_T) /
            2.)
        #_coolant_material_        = flibe.mcf(24)
        _outlet_coolant_material_ = flibe.mcf(25, tmp=self.coolant_outlet_T)
        _Heat_Exchangers_ = self.mat_HX.mcf(
            700,
            tmp=(
                self.coolant_inlet_T +
                self.coolant_outlet_T) /
            2.)
        _control_rod_material_ = control.mcf(50, tmp=self.coolant_inlet_T)
        _control_blade_material_ = control.mcf(
            51,
            tmp=(
                self.coolant_inlet_T +
                self.coolant_outlet_T) /
            2.)
        _shield_absorber_material_ = self.shield.mat_absorber.mcf(
            60,
            tmp=self.coolant_inlet_T)
        _shield_clad_material_ = self.shield.mat_clad.mcf(
            61,
            tmp=self.coolant_inlet_T)

        # Inner Graphite Reflector
        theta = 2*math.pi/self.control_rods
        dt = self.r_control_rods*self.dt_control_rods
        t = self.r_control_rods*self.t_control_rods
        R = self.r_inner_solid_reflector - 1.5*self.r_control_rods
        _inner_reflector_cell_ = '900  1 -_inner_solid_reflector_density_  -901     -602 202'
        _inner_reflector_cell1_ = '901  1 -_inner_solid_reflector_density_  -903 901 -602 402'
        _inner_porous_reflector_cell_ = '902  5 -_density_outer_porous_reflector_ -903 901 -402 202'
        _control_rod_surfaces_ = 'c\nc control rod surfaces \nc'
        _control_rod_cells_ = 'c\nc control rod cells\nc'

        if 2.*math.pi*R/self.control_rods < 3.*self.r_control_rods:
            print('too many control blades or channel is too large might cause overlap or thin webbing between them!!')

        for i in range(self.control_rods):
            x = math.cos(i*theta)*R
            y = math.sin(i*theta)*R

            # generate surfaces
            _control_rod_surfaces_ += '\n%d c/z %1.5E %1.5E %1.5E' % (
                int(8000+i*10), x,  y, self.r_control_rods)
            _inner_reflector_cell_ += '\n          %d' % (int(8000+i*10))
            _inner_reflector_cell1_ += '\n          %d' % (int(8000+i*10))
            _inner_porous_reflector_cell_ += '\n          %d' % (int(8000+i*10))
            # surface 1
            _control_rod_surfaces_ += '\n%d py %1.5E' % (
                int(8000 + i * 10 + 1), (y + t))
            # surface 2
            _control_rod_surfaces_ += '\n%d px %1.5E' % (
                int(8000 + i * 10 + 2), (x + dt))
            # surface 3
            _control_rod_surfaces_ += '\n%d py %1.5E' % (
                int(8000 + i * 10 + 3), (y - t))
            # surface 4
            _control_rod_surfaces_ += '\n%d px %1.5E' % (
                int(8000 + i * 10 + 4), (x - dt))
            # surface 5
            _control_rod_surfaces_ += '\n%d py %1.5E' % (
                int(8000 + i * 10 + 5), (y + dt))
            # surface 6
            _control_rod_surfaces_ += '\n%d px %1.5E' % (
                int(8000 + i * 10 + 6), (x + t))
            # surface 7
            _control_rod_surfaces_ += '\n%d py %1.5E' % (
                int(8000 + i * 10 + 7), (y - dt))
            # surface 8
            _control_rod_surfaces_ += '\n%d px %1.5E' % (
                int(8000 + i * 10 + 8), (x - t))
            _control_rod_cells_ += '\nc CR withdrawn %d 23 -_coolant_inlet_density_ -%d -600 202 imp:n=1 tmp=_coolant_inlet_temperature_' % (
                int(8000+i*10), int(8000+i*10))
            _control_rod_cells_ += '\nc CR withdrawn %d 50 -_control_rod_density_\nc CR withdrawn            ' % (
                int(8000 + i*10 + 1))
            _control_rod_cells_ += '(-%d -%d %d %d -602 600):(-%d -%d %d %d -602 600)' % (int(8000 + i*10 + 1),
                                                                                          int(8000 + i*10 + 2),
                                                                                          int(8000 + i*10 + 3),
                                                                                          int(8000 + i*10 + 4),
                                                                                          int(8000 + i*10 + 5),
                                                                                          int(8000 + i*10 + 6),
                                                                                          int(8000 + i*10 + 7),
                                                                                          int(8000 + i*10 + 8),
                                                                                          )
            _control_rod_cells_ += '\nc CR withdrawn           imp:n=1 tmp=_coolant_inlet_temperature_'
            _control_rod_cells_ += '\nc CR withdrawn %d 23 -_coolant_inlet_density_ -%d #%d -602 600 imp:n=1 tmp=_coolant_inlet_temperature_' % (
                int(8000+i*10+2), int(8000+i*10), int(8000+i*10+1))
            _control_rod_cells_ += '\nc CR engaged %d 23 -_coolant_inlet_density_ -%d #%d -408 406 imp:n=1 tmp=_coolant_inlet_temperature_' % (
                int(8000+i*10), int(8000+i*10), int(8000+i*10+1))
            _control_rod_cells_ += '\nc CR engaged %d 50 -_control_rod_density_\nc CR engaged           ' % (
                int(8000 + i*10 + 1))
            _control_rod_cells_ += '(-%d -%d %d %d -408 406):(-%d -%d %d %d -408 406)' % (int(8000 + i*10 + 1),
                                                                                          int(8000 + i*10 + 2),
                                                                                          int(8000 + i*10 + 3),
                                                                                          int(8000 + i*10 + 4),
                                                                                          int(8000 + i*10 + 5),
                                                                                          int(8000 + i*10 + 6),
                                                                                          int(8000 + i*10 + 7),
                                                                                          int(8000 + i*10 + 8),
                                                                                          )
            _control_rod_cells_ += '\nc CR engaged           imp:n=1 tmp=_coolant_inlet_temperature_'
            _control_rod_cells_ += '\nc CR engaged %d 23 -_coolant_inlet_density_ -%d -602 408 imp:n=1 tmp=_coolant_inlet_temperature_' % (
                int(8000+i*10+3), int(8000+i*10),)
            _control_rod_cells_ += '\nc CR engaged %d 23 -_coolant_inlet_density_ -%d -406 202 imp:n=1 tmp=_coolant_inlet_temperature_' % (
                int(8000+i*10+4), int(8000+i*10),)

        _inner_reflector_cell_ += '\n           imp:n=1 tmp=_coolant_inlet_temperature_\n'
        _inner_reflector_cell1_ += '\n           imp:n=1 tmp=_coolant_inlet_temperature_\n'
        _inner_porous_reflector_cell_ += '\n           imp:n=1 tmp=_coolant_inlet_temperature_'
        _inner_reflector_cell_ += _inner_reflector_cell1_ + \
            _inner_porous_reflector_cell_

        # generate surfaces and cells for control rod channels

        theta = 2*math.pi/self.control_blades
        t = self.dt_control_blades*self.fuel.r_shell*2
        M = self.dm_control_blades*self.fuel.r_shell*2
        N = self.dn_control_blades*self.fuel.r_shell*2
        R = self.r_defueling_pebble_reflector + self.fuel.r_shell * \
            2*self.dm_control_blades + self.dr_control_blades

        _heat_exchanger_reflector_cell_ = '710 6  -_density_outer_porous_reflector_ -701  609 -602 700'
        _defueling_reflector_cell_ = '610 6  -_density_outer_porous_reflector_ -411  609 -700 600'
        _converging_reflector_cell_ = '510 6  -_density_outer_porous_reflector_ -411  509 -600 500'
        _control_blade_surface_ = 'c\nc control blade surfaces \nc'
        _control_blade_cell_ = 'c\nc control blade cells \nc'
        _control_blade_cell2_ = ''
        _CB_cells_ = ''

        for i in range(self.control_blades):
            x = math.cos(i*theta)*R
            y = math.sin(i*theta)*R

            m = [math.cos(i*theta), math.sin(i*theta)]
            n = [math.cos(i*theta+math.pi/2.), math.sin(i*theta+math.pi/2.)]

            # point A
            # north east corner on n-direction blade
            #

            A = [(R*m[0] + t*m[0] + N*n[0]), (R*m[1] + t*m[1] + N*n[1])]

            # point B
            # south west corner on n-direction blade

            B = [(R*m[0] - t*m[0] - N*n[0]), (R*m[1] - t*m[1] - N*n[1])]

            # point C
            # north east corner on m-direction blade

            C = [(R*m[0] + t*n[0] + M*m[0]), (R*m[1] + t*n[1] + M*m[1])]

            # point D
            # north east corder on m direction blade

            D = [(R*m[0] - t*n[0] - M*m[0]), (R*m[1] - t*n[1] - M*m[1])]

            _control_blade_cell_ += '\nc CB withdrawn %d 25 -_coolant_outlet_density_ \nc CB withdrawn          (' % (
                9000 + i*10)
            _control_blade_cell2_ += '\nc CB engaged %d 51 -_control_blade_density_ \nc CB engaged          (' % (
                9000 + i*10)
            _defueling_reflector_cell_ += '\n          ('
            _converging_reflector_cell_ += '\n          ('
            _heat_exchanger_reflector_cell_ += '\n          ('

            # surface 1

            b = A[0]*math.tan(i*theta) - A[1]

            if (x*math.tan(i*theta) - 1*y - b) < 0:
                # sense of center [x,y] is negative with respect to the surface
                _control_blade_cell_ += '-%d ' % (9000 + i*10 + 1)
                _control_blade_cell2_ += '-%d ' % (9000 + i*10 + 1)
                _converging_reflector_cell_ += '%d:' % (9000 + i*10 + 1)
                _defueling_reflector_cell_ += '%d:' % (9000 + i*10 + 1)
                _heat_exchanger_reflector_cell_ += '%d:' % (9000 + i*10 + 1)
            else:
                _control_blade_cell_ += '%d ' % (9000 + i*10 + 1)
                _control_blade_cell2_ += '%d ' % (9000 + i*10 + 1)
                _converging_reflector_cell_ += '-%d:' % (9000 + i*10 + 1)
                _defueling_reflector_cell_ += '-%d:' % (9000 + i*10 + 1)
                _heat_exchanger_reflector_cell_ += '-%d:' % (9000 + i*10 + 1)

            _control_blade_surface_ += '\n%d p %1.5E -1 0 %1.5E' % (
                (9000 + i*10 + 1), math.tan(i*theta), b)

            # surface 2

            b = A[0]*math.tan(i*theta + math.pi/2.) - A[1]

            if (x*math.tan(i*theta + math.pi/2) - 1*y - b) < 0:
                # sense of center [x,y] is negative with respect to the surface
                _control_blade_cell_ += '-%d ' % (9000 + i*10 + 2)
                _control_blade_cell2_ += '-%d ' % (9000 + i*10 + 2)
                _converging_reflector_cell_ += '%d:' % (9000 + i*10 + 2)
                _defueling_reflector_cell_ += '%d:' % (9000 + i*10 + 2)
                _heat_exchanger_reflector_cell_ += '%d:' % (9000 + i*10 + 2)
            else:
                _control_blade_cell_ += '%d ' % (9000 + i*10 + 2)
                _control_blade_cell2_ += '%d ' % (9000 + i*10 + 2)
                _converging_reflector_cell_ += '-%d:' % (9000 + i*10 + 2)
                _defueling_reflector_cell_ += '-%d:' % (9000 + i*10 + 2)
                _heat_exchanger_reflector_cell_ += '-%d:' % (9000 + i*10 + 2)

            _control_blade_surface_ += '\n%d p %1.5E -1 0 %1.5E' % (
                (9000 + i*10 + 2), math.tan(i*theta + math.pi/2.), b)

            # surface 3

            b = B[0]*math.tan(i*theta) - B[1]

            if (x*math.tan(i*theta) - 1*y - b) < 0:
                # sense of center [x,y] is negative with respect to the surface
                _control_blade_cell_ += '-%d ' % (9000 + i*10 + 3)
                _control_blade_cell2_ += '-%d ' % (9000 + i*10 + 3)
                _converging_reflector_cell_ += '%d:' % (9000 + i*10 + 3)
                _defueling_reflector_cell_ += '%d:' % (9000 + i*10 + 3)
                _heat_exchanger_reflector_cell_ += '%d:' % (9000 + i*10 + 3)
            else:
                _control_blade_cell_ += '%d ' % (9000 + i*10 + 3)
                _control_blade_cell2_ += '%d ' % (9000 + i*10 + 3)
                _converging_reflector_cell_ += '-%d:' % (9000 + i*10 + 3)
                _defueling_reflector_cell_ += '-%d:' % (9000 + i*10 + 3)
                _heat_exchanger_reflector_cell_ += '-%d:' % (9000 + i*10 + 3)

            _control_blade_surface_ += '\n%d p %1.5E -1 0 %1.5E' % (
                (9000 + i*10 + 3), math.tan(i*theta), b)

            # surface 4

            b = B[0]*math.tan(i*theta + math.pi/2.) - B[1]

            if (x*math.tan(i*theta + math.pi/2.) - 1*y - b) < 0:
                # sense of center [x,y] is negative with respect to the surface
                _control_blade_cell_ += '-%d ' % (9000 + i*10 + 4)
                _control_blade_cell2_ += '-%d ' % (9000 + i*10 + 4)
                _converging_reflector_cell_ += '%d:' % (9000 + i*10 + 4)
                _defueling_reflector_cell_ += '%d:' % (9000 + i*10 + 4)
                _heat_exchanger_reflector_cell_ += '%d:' % (9000 + i*10 + 4)
            else:
                _control_blade_cell_ += '%d ' % (9000 + i*10 + 4)
                _control_blade_cell2_ += '%d ' % (9000 + i*10 + 4)
                _converging_reflector_cell_ += '-%d:' % (9000 + i*10 + 4)
                _defueling_reflector_cell_ += '-%d:' % (9000 + i*10 + 4)
                _heat_exchanger_reflector_cell_ += '-%d:' % (9000 + i*10 + 4)

            _control_blade_surface_ += '\n%d p %1.5E -1 0 %1.5E' % (
                (9000 + i*10 + 4), math.tan(i*theta + math.pi/2.), b)

            # generated all the surfaces for the wing in the N direction!

            _control_blade_cell_ = _control_blade_cell_[:-1] + ' -602 509):('
            _control_blade_cell2_ = _control_blade_cell2_[:-1] + ' -602 400):('
            _converging_reflector_cell_ = _converging_reflector_cell_[
                : -1] + ') ('
            _defueling_reflector_cell_ = _defueling_reflector_cell_[:-1] + ') ('
            _heat_exchanger_reflector_cell_ = _heat_exchanger_reflector_cell_[
                : -1] + ') ('

            # surface 5

            b = C[0]*math.tan(i*theta) - C[1]

            if (x*math.tan(i*theta) - 1*y - b) < 0:
                # sense of center [x,y] is negative with respect to the surface
                _control_blade_cell_ += '-%d ' % (9000 + i*10 + 5)
                _control_blade_cell2_ += '-%d ' % (9000 + i*10 + 5)
                _converging_reflector_cell_ += '%d:' % (9000 + i*10 + 5)
                _defueling_reflector_cell_ += '%d:' % (9000 + i*10 + 5)
                _heat_exchanger_reflector_cell_ += '%d:' % (9000 + i*10 + 5)
            else:
                _control_blade_cell_ += '%d ' % (9000 + i*10 + 5)
                _control_blade_cell2_ += '%d ' % (9000 + i*10 + 5)
                _converging_reflector_cell_ += '-%d:' % (9000 + i*10 + 5)
                _defueling_reflector_cell_ += '-%d:' % (9000 + i*10 + 5)
                _heat_exchanger_reflector_cell_ += '-%d:' % (9000 + i*10 + 5)

            _control_blade_surface_ += '\n%d p %1.5E -1 0 %1.5E' % (
                (9000 + i*10 + 5), math.tan(i*theta), b)

            # surface 6

            b = C[0]*math.tan(i*theta + math.pi/2.) - C[1]

            if (x*math.tan(i*theta + math.pi/2.) - 1*y - b) < 0:
                # sense of center [x,y] is negative with respect to the surface
                _control_blade_cell_ += '-%d ' % (9000 + i*10 + 6)
                _control_blade_cell2_ += '-%d ' % (9000 + i*10 + 6)
                _converging_reflector_cell_ += '%d:' % (9000 + i*10 + 6)
                _defueling_reflector_cell_ += '%d:' % (9000 + i*10 + 6)
                _heat_exchanger_reflector_cell_ += '%d:' % (9000 + i*10 + 6)
            else:
                _control_blade_cell_ += '%d ' % (9000 + i*10 + 6)
                _control_blade_cell2_ += '%d ' % (9000 + i*10 + 6)
                _converging_reflector_cell_ += '-%d:' % (9000 + i*10 + 6)
                _defueling_reflector_cell_ += '-%d:' % (9000 + i*10 + 6)
                _heat_exchanger_reflector_cell_ += '-%d:' % (9000 + i*10 + 6)

            _control_blade_surface_ += '\n%d p %1.5E -1 0 %1.5E' % (
                (9000 + i*10 + 6), math.tan(i*theta + math.pi/2.), b)

            # surface 7

            b = D[0]*math.tan(i*theta) - D[1]

            if (x*math.tan(i*theta) - 1*y - b) < 0:
                # sense of center [x,y] is negative with respect to the surface
                _control_blade_cell_ += '-%d ' % (9000 + i*10 + 7)
                _control_blade_cell2_ += '-%d ' % (9000 + i*10 + 7)
                _converging_reflector_cell_ += '%d:' % (9000 + i*10 + 7)
                _defueling_reflector_cell_ += '%d:' % (9000 + i*10 + 7)
                _heat_exchanger_reflector_cell_ += '%d:' % (9000 + i*10 + 7)
            else:
                _control_blade_cell_ += '%d ' % (9000 + i*10 + 7)
                _control_blade_cell2_ += '%d ' % (9000 + i*10 + 7)
                _converging_reflector_cell_ += '-%d:' % (9000 + i*10 + 7)
                _defueling_reflector_cell_ += '-%d:' % (9000 + i*10 + 7)
                _heat_exchanger_reflector_cell_ += '-%d:' % (9000 + i*10 + 7)

            _control_blade_surface_ += '\n%d p %1.5E -1 0 %1.5E' % (
                (9000 + i*10 + 7), math.tan(i*theta), b)

            # surface 8

            b = D[0]*math.tan(i*theta + math.pi/2.) - D[1]

            if (x*math.tan(i*theta + math.pi/2.) - 1*y - b) < 0:
                # sense of center [x,y] is negative with respect to the surface
                _control_blade_cell_ += '-%d ' % (9000 + i*10 + 8)
                _control_blade_cell2_ += '-%d ' % (9000 + i*10 + 8)
                _converging_reflector_cell_ += '%d:' % (9000 + i*10 + 8)
                _defueling_reflector_cell_ += '%d:' % (9000 + i*10 + 8)
                _heat_exchanger_reflector_cell_ += '%d:' % (9000 + i*10 + 8)
            else:
                _control_blade_cell_ += '%d ' % (9000 + i*10 + 8)
                _control_blade_cell2_ += '%d ' % (9000 + i*10 + 8)
                _converging_reflector_cell_ += '-%d:' % (9000 + i*10 + 8)
                _defueling_reflector_cell_ += '-%d:' % (9000 + i*10 + 8)
                _heat_exchanger_reflector_cell_ += '-%d:' % (9000 + i*10 + 8)

            _control_blade_surface_ += '\n%d p %1.5E -1 0 %1.5E' % (
                (9000 + i*10 + 8), math.tan(i*theta + math.pi/2.), b)

            # generated all the surfaces for the wing in the N direction!

            _control_blade_cell_ = _control_blade_cell_[
                :-1] + ' -602 509)\nc CB withdrawn          imp:n=1 tmp=_coolant_outlet_temperature_'
            _control_blade_cell2_ = _control_blade_cell2_[
                :-1] + ' -602 400)\nc CB engaged          imp:n=1 tmp=_coolant_temperature_'
            _converging_reflector_cell_ = _converging_reflector_cell_[:-1] + ')'
            _defueling_reflector_cell_ = _defueling_reflector_cell_[:-1] + ')'
            _heat_exchanger_reflector_cell_ = _heat_exchanger_reflector_cell_[
                : -1] + ')'
            _CB_cells_ += '\nc CB engaged           #%d' % (int(9000+i*10))

        _converging_reflector_cell_ += '\n          imp:n=1 tmp=_coolant_outlet_temperature_'
        _defueling_reflector_cell_ += '\n          imp:n=1 tmp=_coolant_outlet_temperature_'
        _heat_exchanger_reflector_cell_ += '\n          imp:n=1 tmp=_coolant_outlet_temperature_'

        _control_rod_density_ = '%1.5E' % control.mass()
        _control_blade_density_ = '%1.5E' % control.mass()

        # generate cells and surfaces for control blade channels

        # Replace tokens in skeleton file
        spec_tally = '\nfc4 Neutron Spectrum Tally\nf4:n      ' + _fuel_cells_[
            10:] + '         (' + _fuel_cells_[
            10:-1] + ')\n e4  1.0000E-11 1.0000E-10 5.0000E-10 7.5000E-10 1.0000E-09 1.2000E-09\n     1.5000E-09 2.0000E-09 2.5000E-09 3.0000E-09 4.0000E-09 5.0000E-09\n     7.5000E-09 1.0000E-08 2.5300E-08 3.0000E-08 4.0000E-08 5.0000E-08\n     6.0000E-08 7.0000E-08 8.0000E-08 9.0000E-08 1.0000E-07 1.2500E-07\n     1.5000E-07 1.7500E-07 2.0000E-07 2.2500E-07 2.5000E-07 2.7500E-07\n     3.0000E-07 3.2500E-07 3.5000E-07 3.7500E-07 4.0000E-07 4.5000E-07\n     5.0000E-07 5.5000E-07 6.0000E-07 6.2500E-07 6.5000E-07 7.0000E-07\n     7.5000E-07 8.0000E-07 8.5000E-07 9.0000E-07 9.2500E-07 9.5000E-07\n     9.7500E-07 1.0000E-06 1.0100E-06 1.0200E-06 1.0300E-06 1.0400E-06\n     1.0500E-06 1.0600E-06 1.0700E-06 1.0800E-06 1.0900E-06 1.1000E-06\n     1.1100E-06 1.1200E-06 1.1300E-06 1.1400E-06 1.1500E-06 1.1750E-06\n     1.2000E-06 1.2250E-06 1.2500E-06 1.3000E-06 1.3500E-06 1.4000E-06\n     1.4500E-06 1.5000E-06 1.5900E-06 1.6800E-06 1.7700E-06 1.8600E-06\n     1.9400E-06 2.0000E-06 2.1200E-06 2.2100E-06 2.3000E-06 2.3800E-06\n     2.4700E-06 2.5700E-06 2.6700E-06 2.7700E-06 2.8700E-06 2.9700E-06\n     3.0000E-06 3.0500E-06 3.1500E-06 3.5000E-06 3.7300E-06 4.0000E-06\n     4.7500E-06 5.0000E-06 5.4000E-06 6.0000E-06 6.2500E-06 6.5000E-06\n     6.7500E-06 7.0000E-06 7.1500E-06 8.1000E-06 9.1000E-06 1.0000E-05\n     1.1500E-05 1.1900E-05 1.2900E-05 1.3750E-05 1.4400E-05 1.5100E-05\n     1.6000E-05 1.7000E-05 1.8500E-05 1.9000E-05 2.0000E-05 2.1000E-05\n     2.2500E-05 2.5000E-05 2.7500E-05 3.0000E-05 3.1250E-05 3.1750E-05\n     3.3250E-05 3.3750E-05 3.4600E-05 3.5500E-05 3.7000E-05 3.8000E-05\n     3.9100E-05 3.9600E-05 4.1000E-05 4.2400E-05 4.4000E-05 4.5200E-05\n     4.7000E-05 4.8300E-05 4.9200E-05 5.0600E-05 5.2000E-05 5.3400E-05\n     5.9000E-05 6.1000E-05 6.5000E-05 6.7500E-05 7.2000E-05 7.6000E-05\n     8.0000E-05 8.2000E-05 9.0000E-05 1.0000E-04 1.0800E-04 1.1500E-04\n     1.1900E-04 1.2200E-04 1.8600E-04 1.9250E-04 2.0750E-04 2.1000E-04\n     2.4000E-04 2.8500E-04 3.0500E-04 5.5000E-04 6.7000E-04 6.8300E-04\n     9.5000E-04 1.1500E-03 1.5000E-03 1.5500E-03 1.8000E-03 2.2000E-03\n     2.2900E-03 2.5800E-03 3.0000E-03 3.7400E-03 3.9000E-03 6.0000E-03\n     8.0300E-03 9.5000E-03 1.3000E-02 1.7000E-02 2.5000E-02 3.0000E-02\n     4.5000E-02 5.0000E-02 5.2000E-02 6.0000E-02 7.3000E-02 7.5000E-02\n     8.2000E-02 8.5000E-02 1.0000E-01 1.2830E-01 1.5000E-01 2.0000E-01\n     2.7000E-01 3.3000E-01 4.0000E-01 4.2000E-01 4.4000E-01 4.7000E-01\n     4.9952E-01 5.5000E-01 5.7300E-01 6.0000E-01 6.7000E-01 6.7900E-01\n     7.5000E-01 8.2000E-01 8.6110E-01 8.7500E-01 9.0000E-01 9.2000E-01\n     1.0100E+00 1.1000E+00 1.2000E+00 1.2500E+00 1.3170E+00 1.3560E+00\n     1.4000E+00 1.5000E+00 1.8500E+00 2.3540E+00 2.4790E+00 3.0000E+00\n     4.3040E+00 4.8000E+00 6.4340E+00 8.1873E+00 1.0000E+01 1.2840E+01\n     1.3840E+01 1.4550E+01 1.5683E+01 1.7333E+01 2.0000E+01\nsd4 1 ' + str(
            int(
                4*5*self.passes)) + 'r'
        dose_tally = '\n'

        # tally of dose to inner graphite reflector
        dose_tally += 'fc114 fast flux mesh for central graphite reflector'
        dose_tally += '\nfmesh114:n geom=cyl'
        dose_tally += '\n          imesh=%1.5E %1.5E iints=1 10' % (
            (self.r_inner_solid_reflector - self.inner_porous_reflector_dr),
            self.r_inner_solid_reflector)
        dose_tally += '\n          jmesh=%1.5E %1.5E jints=1 20' % (
            (self.dz_entrance + self.dz_expansion),
            (self.dz_entrance + self.dz_expansion + self.dz_active))
        dose_tally += '\n          kmesh=1                       kints=1'
        dose_tally += '\n          emesh=     .1          20'

        if self.zf_inner_porous > ((self.dz_expansion + self.dz_active)/(self.dz_expansion + self.dz_active + self.dz_converging)):
            # the porous region is continuous in the active region
            dose_tally += '\nfc124 DPA mesh for central graphite reflector (porous)'
            dose_tally += '\nfmesh124:n geom=cyl'
            dose_tally += '\n          imesh=%1.5E %1.5E iints=1 10' % (
                (self.r_inner_solid_reflector - self.inner_porous_reflector_dr),
                self.r_inner_solid_reflector)
            dose_tally += '\n          jmesh=%1.5E %1.5E jints=1 20' % (
                (self.dz_entrance + self.dz_expansion),
                (self.dz_entrance + self.dz_expansion + self.dz_active))
            dose_tally += '\n          kmesh=1                       kints=1'
            dose_tally += '\nfm124     (1 5 444)'

        else:
            # there is a porous region / solid region interface in the active
            # region
            dose_tally += '\nfc124 fast flux mesh for central graphite reflector (porous)'
            dose_tally += '\nfmesh124:n geom=cyl'
            dose_tally += '\n          imesh=%1.5E %1.5E iints=1 10' % (
                (self.r_inner_solid_reflector - self.inner_porous_reflector_dr),
                self.r_inner_solid_reflector)
            dose_tally += '\n          jmesh=%1.5E %1.5E jints=1 20' % ((self.dz_entrance + self.dz_expansion), (
                self.dz_entrance + (self.dz_expansion + self.dz_active + self.dz_converging)*self.zf_inner_porous))
            dose_tally += '\n          kmesh=1                       kints=1'
            dose_tally += '\nfm124     (1 5 444)'

            dose_tally += '\nfc134 DPA mesh for central graphite reflector (solid)'
            dose_tally += '\nfmesh134:n geom=cyl'
            dose_tally += '\n          imesh=%1.5E %1.5E iints=1 10' % (
                (self.r_inner_solid_reflector - self.inner_porous_reflector_dr),
                self.r_inner_solid_reflector)
            dose_tally += '\n          jmesh=%1.5E %1.5E jints=1 20' % ((self.dz_entrance + (
                self.dz_expansion + self.dz_active + self.dz_converging)*self.zf_inner_porous),
                (self.dz_entrance + self.dz_expansion + self.dz_active))
            dose_tally += '\n          kmesh=1                       kints=1'
            dose_tally += '\nfm134     (1 1 444)'

        # tally of dose to outer graphite reflector
        dose_tally += '\nfc144 fast flux mesh for outer solid graphite reflector'
        dose_tally += '\nfmesh144:n geom=cyl'
        dose_tally += '\n          imesh=%1.5E %1.5E iints=1 10' % (
            self.r_active_pebble_reflector,
            (self.r_active_pebble_reflector + self.outer_porous_reflector_dr))
        dose_tally += '\n          jmesh=%1.5E %1.5E jints=1 20' % (
            (self.dz_entrance + self.dz_expansion),
            (self.dz_entrance + self.dz_expansion + self.dz_active))
        dose_tally += '\n          kmesh=1                       kints=1'
        dose_tally += '\n          emesh=     .1          20'

        dose_tally += '\nfc154 DPA mesh for outer solid graphite reflector (solid)'
        dose_tally += '\nfmesh154:n geom=cyl'
        dose_tally += '\n          imesh=%1.5E %1.5E iints=1 10' % (
            self.r_active_pebble_reflector,
            (self.r_active_pebble_reflector + self.outer_porous_reflector_dr))
        dose_tally += '\n          jmesh=%1.5E %1.5E jints=1 20' % (
            (self.dz_entrance + self.dz_expansion),
            (self.dz_entrance + self.dz_expansion + self.dz_active
             * (1. - .999 * self.zf_inner_porous)))
        dose_tally += '\n          kmesh=1                       kints=1'
        dose_tally += '\nfm154     (1 2 444)'

        dose_tally += '\nfc164 DPA mesh for outer solid graphite reflector (porous)'
        dose_tally += '\nfmesh164:n geom=cyl'
        dose_tally += '\n          imesh=%1.5E %1.5E iints=1 10' % (
            self.r_active_pebble_reflector,
            (self.r_active_pebble_reflector + self.outer_porous_reflector_dr))
        dose_tally += '\n          jmesh=%1.5E %1.5E jints=1 20' % (
            (self.dz_entrance + self.dz_expansion + self.dz_active
             * (1. - .999 * self.zf_inner_porous)),
            (self.dz_entrance + self.dz_expansion + self.dz_active))
        dose_tally += '\n          kmesh=1                       kints=1'
        dose_tally += '\nfm164     (1 6 444)'

        # tally of dose to heat exchanger
        dose_tally += '\nfc174 fast flux mesh for heat exchanger'
        dose_tally += '\nfmesh174:n geom=cyl'
        dose_tally += '\n          imesh=%1.5E %1.5E iints=1 10' % (
            (
                (self.r_defueling_pebble_reflector + self.fuel.r_shell * 2 * 2 *
                 self.dm_control_blades + 2 * self.dr_control_blades)),
            self.r_core_barrel)
        dose_tally += '\n          jmesh=%1.5E %1.5E jints=1 10' % (
            (self.dz_entrance + self.dz_expansion + self.dz_active +
             self.dz_converging + self.dz_HX_clearance), self.Z)
        dose_tally += '\n          kmesh=1                       kints=1'
        dose_tally += '\n          emesh=     .1          20'

        dose_tally += '\nfc184 DPA mesh for heat exchanger'
        dose_tally += '\nfmesh184:n geom=cyl'
        dose_tally += '\n          imesh=%1.5E %1.5E iints=1 10' % (
            (
                (self.r_defueling_pebble_reflector + self.fuel.r_shell * 2 * 2 *
                 self.dm_control_blades + 2 * self.dr_control_blades)),
            self.r_core_barrel)
        dose_tally += '\n          jmesh=%1.5E %1.5E jints=1 10' % (
            (self.dz_entrance + self.dz_expansion + self.dz_active +
             self.dz_converging + self.dz_HX_clearance), self.Z)
        dose_tally += '\n          kmesh=1                       kints=1'
        dose_tally += '\nfm184 (1 700 444)'

        # tally of dose to core barrel
        dose_tally += '\nfc194 fast flux mesh for core barrel'
        dose_tally += '\nfmesh194:n geom=cyl'
        dose_tally += '\n          imesh=%1.5E %1.5E iints=1 10' % (
            self.r_active_solid_reflector, self.r_core_barrel)
        dose_tally += '\n          jmesh=%1.5E %1.5E jints=1 20' % (
            (self.dz_entrance + self.dz_expansion),
            (self.dz_entrance + self.dz_expansion + self.dz_active))
        dose_tally += '\n          kmesh=1                       kints=1'
        dose_tally += '\n          emesh=     .1          20'

        dose_tally += '\nfc204 DPA mesh for core barrel'
        dose_tally += '\nfmesh204:n geom=cyl'
        dose_tally += '\n          imesh=%1.5E %1.5E iints=1 10' % (
            self.r_active_solid_reflector, self.r_core_barrel)
        dose_tally += '\n          jmesh=%1.5E %1.5E jints=1 20' % (
            (self.dz_entrance + self.dz_expansion),
            (self.dz_entrance + self.dz_expansion + self.dz_active))
        dose_tally += '\n          kmesh=1                       kints=1'
        dose_tally += '\nfm204 (1 1 444)'

        # tally of dose to reactor vessel
        dose_tally += '\nfc214 fast flux mesh for reactor vessel'
        dose_tally += '\nfmesh214:n geom=cyl'
        dose_tally += '\n          imesh=%1.5E %1.5E iints=1 1' % (
            self.r_HX, self.r_reactor_vessel)
        dose_tally += '\n          jmesh=%1.5E %1.5E jints=1 20' % (
            (self.dz_entrance + self.dz_expansion),
            (self.dz_entrance + self.dz_expansion + self.dz_active))
        dose_tally += '\n          kmesh=1                       kints=1'
        dose_tally += '\n          emesh=     .1          20'

        dose_tally += '\nfc224 DPA mesh for reactor vessel'
        dose_tally += '\nfmesh224:n geom=cyl'
        dose_tally += '\n          imesh=%1.5E %1.5E iints=1 1' % (
            self.r_HX, self.r_reactor_vessel)
        dose_tally += '\n          jmesh=%1.5E %1.5E jints=1 20' % (
            (self.dz_entrance + self.dz_expansion),
            (self.dz_entrance + self.dz_expansion + self.dz_active))
        dose_tally += '\n          kmesh=1                       kints=1'
        dose_tally += '\nfm224 (1 31 444)'

        # tally of dose to pebble separator
        dose_tally += '\nfc234 fast flux mesh for pebble separator'
        dose_tally += '\nfmesh234:n geom=cyl'
        dose_tally += '\n          imesh=%1.5E %1.5E iints=1 1' % (
            self.r_entrance_fuel-.5, self.r_entrance_fuel+.5)
        dose_tally += '\n          jmesh=%1.5E %1.5E jints=1 20' % (
            1E-30, (self.dz_entrance))
        dose_tally += '\n          kmesh=1                       kints=1'
        dose_tally += '\n          emesh=     .1          20'

        dose_tally += '\nfc244 dummy tally to turn on fm capability\nf244:n 900\nfm244 (1 1 444)\nsd244 1'

        # economy tallies

        # low burnup kernel

        economy_tally = '\n'

        economy_tally += 'fc304 fuel kernel tallies\n'
        economy_tally += ' f304:n   ' + _fuel_cells_[10:-1] + '\n'
        economy_tally += 'fm304:n\n' + _fuel_tally_
        economy_tally += 'sd304 1 %dr\n' % (4*5*self.passes - 1)

        economy_tally += 'fc314 TRISO layers + matrix\n'
        economy_tally += ' f314:n  (' + _matrix_cells_[10:-1] + ')\n'
        economy_tally += 'fm314 (-1 11 (-2))\nsd314 1\n'

        economy_tally += 'fc324 pebble cores\n'
        economy_tally += ' f324:n  (' + _pebble_core_cells_[10:-1] + ')\n'
        economy_tally += 'fm324 (-1 3 (-2))\nsd324 1\n'

        economy_tally += 'fc334 pebble shells\n'
        economy_tally += ' f334:n  (' + _pebble_shell_cells_[10:-1] + ')\n'
        economy_tally += 'fm334 (-1 4 (-2))\nsd334 1\n'

        economy_tally += 'fc344 coolant\n'
        economy_tally += ' f344:n  (' + _coolant_cells_[10:-1] + ')\n'
        economy_tally += 'fm344 (-1 24 (-2) (16) (17))\nsd344 1'

        # leakage tallies

        # net leakage in radially to inner reflector

        leakage_tally = '\nc\nfc901 net leakage to inner reflector'
        leakage_tally += '\n f901:n 903'
        leakage_tally += '\n c901 0 1'
        leakage_tally += '\nfs901 -204 602'
        leakage_tally += '\nsd901 1 1 1\nc'

        # net leakage out axially to lower plenum

        leakage_tally += '\nfc911 net leakage out axialy to lower plenum'
        leakage_tally += '\n f911:n 204'
        leakage_tally += '\n c911 0 1'
        leakage_tally += '\nfs911 -903 207'
        leakage_tally += '\nsd911 1 1 1\nc'

        # net leakage out radially of entrance region

        leakage_tally += '\nfc921 net leakage out radially of entrance region'
        leakage_tally += '\n f921:n 207'
        leakage_tally += '\n c921 0 1'
        leakage_tally += '\nfs921 -204 300'
        leakage_tally += '\nsd921 1 1 1\nc'

        # net leakage out radiall of expansion region

        leakage_tally += '\nfc931 net leakage out radially of expansion region'
        leakage_tally += '\n f931:n 307'
        leakage_tally += '\n c931 0 1'
        leakage_tally += '\nfs931 400 -300'
        leakage_tally += '\nsd931 1 1 1\nc'

        # net leakage out radially of active region

        leakage_tally += '\nfc941 net leakage out radially of active region'
        leakage_tally += '\n f941:n 407'
        leakage_tally += '\n c941 0 1'
        leakage_tally += '\nfs941 500 -400'
        leakage_tally += '\nsd941 1 1 1\nc'

        # net leakage out radially of converging region

        leakage_tally += '\nfc951 net leakage out radially of converging region'
        leakage_tally += '\n f951:n 507'
        leakage_tally += '\n c951 0 1'
        leakage_tally += '\nfs951 600 -500'
        leakage_tally += '\nsd951 1 1 1\nc'

        # net leakage out radially of defueling chute

        leakage_tally += '\nfc961 net leakage out radially of defueling chute'
        leakage_tally += '\n f961:n 607'
        leakage_tally += '\n c961 0 1'
        leakage_tally += '\nfs961 602 -600'
        leakage_tally += '\nsd961 1 1 1\nc'

        # net leakage out axially to vacuum

        leakage_tally += '\nfc971 net leakage out axially to vacuum '
        leakage_tally += '\n f971:n 602'
        leakage_tally += '\n c971 0 1'
        leakage_tally += '\nfs971 -903 607'
        leakage_tally += '\nsd971 1 1 1\nc'

        # tallies

        _tally_ = ''
        for tal in self.tallies:
            if tal[0] in 'sS':
                _tally_ += spec_tally
            elif tal[0] in 'mM':
                fuel.mocup_strings()
                # determined fuel cells
                _tally_ += fuel.tally.replace('_cells_', _fuel_cells_[10:-1])
                _fuel_materials_ += fuel.single_mat + '\n'
                # tallies for AUSCE
                # determined coolant cells
                _tally_ += '\n f64:n   (' + _coolant_cells_[
                    10:-1] + ')\nfm64 (1 536 -2)\nsd64 1\n f74:n   (' + _coolant_cells_[
                    10:-1] + ')\nfm74 (1 549 107)\nsd74 1'
            elif tal[0] in 'dD':
                _tally_ += dose_tally
            elif tal[0] in 'lL':
                _tally_ += leakage_tally
            elif tal[0] in 'eE':
                _tally_ += economy_tally
                # ensure leakage tallies are in the deck
                flag = 0
                for tal1 in self.tallies:
                    if tal1[0] in 'lL':
                        flag = 1
                if flag == 0:
                    _tally_ += leakage_tally

        _power_ = '%1.1f' % self.power
        _power_density_ = '%1.1f' % self.fuel_power_density
        _fuel_fraction_ = '%1.1f' % (100*(self.r_active_fuel**2 - self.r_inner_solid_reflector**2.) /
                                     (self.r_active_pebble_reflector**2. - self.r_inner_solid_reflector**2.))
        _outer_graphite_reflector_thickness_ = '%1.1f' % (
            self.dr_active_solid_reflector)
        _CHM_ = '%1.1f' % (self.fuel.CHM()[0])
        _kernel_diameter_ = '%1.1f' % (2*10000*self.fuel.TRISO.dr_kernel)
        _pebble_density_ = '%1.3f' % (self.fuel.CHM()[1])

        _absorber_radius_ = '%1.5E' % self.shield.absorber_r
        _absorber_clad_radius_ = '%1.5E' % self.shield.clad_r
        _absorber_hpitch_ = '%1.5E' % (self.shield.pitch/2.)
        _absorber_pitch_ = '%1.5E' % self.shield.pitch
        _shield_inner_radius_ = '%1.5E' % (
            self.r_active_solid_reflector - self.shield_dr)
        _absorber_density_ = '%1.5E' % self.shield.mat_absorber.mass()
        _shield_clad_density_ = '%1.5E' % self.shield.mat_clad.mass()
        _shield_block_density_ = '%1.5E' % self.shield.mat_matrix.mass()

        _wall_pf_ = '%2.1f' % (100*self.wall_PF)
        _bulk_pf_ = '%2.1f' % (100*self.bulk_PF)

        skele = open('skeleton2').read()

        if self.passes == 1:
            pebble_cells_loc = 'pebble_cells1'
            _pebble_cells_ = open(pebble_cells_loc).read()
        elif self.passes == 2:
            pebble_cells_loc = 'pebble_cells2'
            _pebble_cells_ = open(pebble_cells_loc).read()
        elif self.passes == 4:
            pebble_cells_loc = 'pebble_cells4'
            _pebble_cells_ = open(pebble_cells_loc).read()
        elif self.passes == 8:
            pebble_cells_loc = 'pebble_cells8'
            _pebble_cells_ = open(pebble_cells_loc).read()
        else:
            print('there is no pebble cells skeleton for the desired number of passes')
            poop

        skele = skele.replace('_pebble_cells_', _pebble_cells_)
        skele = skele.replace('_inner_reflector_cell_', _inner_reflector_cell_)
        skele = skele.replace('_control_rod_cells_', _control_rod_cells_)
        skele = skele.replace('_control_rod_surfaces_', _control_rod_surfaces_)

        skele = skele.replace('_control_blade_cell_', _control_blade_cell_)
        skele = skele.replace('_control_blade_cell2_', _control_blade_cell2_)
        skele = skele.replace(
            '_control_blade_surface_',
            _control_blade_surface_)
        skele = skele.replace(
            '_converging_reflector_cell_',
            _converging_reflector_cell_)
        skele = skele.replace(
            '_defueling_reflector_cell_',
            _defueling_reflector_cell_)
        skele = skele.replace(
            '_heat_exchanger_reflector_cell_',
            _heat_exchanger_reflector_cell_)
        skele = skele.replace('_CB_cells_', _CB_cells_)

        skele = skele.replace('_fuel_power_density_', _power_density_)
        skele = skele.replace('_power_', _power_)
        skele = skele.replace('_fuel_fraction_', _fuel_fraction_)
        skele = skele.replace(
            '_outer_graphite_reflector_thickness_',
            _outer_graphite_reflector_thickness_)
        skele = skele.replace('_CHM_', _CHM_)
        skele = skele.replace('_kernel_diameter_', _kernel_diameter_)
        skele = skele.replace('_pebble_density_', _pebble_density_)

        skele = skele.replace('_fuel_density_', _fuel_density_)
        skele = skele.replace('_fuel_temperature_', _fuel_temperature_)

        skele = skele.replace('_fuel_volume11_', _fuel_volume11_)
        skele = skele.replace('_fuel_volume21_', _fuel_volume21_)
        skele = skele.replace('_fuel_volume31_', _fuel_volume31_)
        skele = skele.replace('_fuel_volume41_', _fuel_volume41_)

        skele = skele.replace('_fuel_volume12_', _fuel_volume12_)
        skele = skele.replace('_fuel_volume22_', _fuel_volume22_)
        skele = skele.replace('_fuel_volume32_', _fuel_volume32_)
        skele = skele.replace('_fuel_volume42_', _fuel_volume42_)

        skele = skele.replace('_fuel_volume13_', _fuel_volume13_)
        skele = skele.replace('_fuel_volume23_', _fuel_volume23_)
        skele = skele.replace('_fuel_volume33_', _fuel_volume33_)
        skele = skele.replace('_fuel_volume43_', _fuel_volume43_)

        skele = skele.replace('_fuel_volume14_', _fuel_volume14_)
        skele = skele.replace('_fuel_volume24_', _fuel_volume24_)
        skele = skele.replace('_fuel_volume34_', _fuel_volume34_)
        skele = skele.replace('_fuel_volume44_', _fuel_volume44_)

        skele = skele.replace('_fuel_volume15_', _fuel_volume15_)
        skele = skele.replace('_fuel_volume25_', _fuel_volume25_)
        skele = skele.replace('_fuel_volume35_', _fuel_volume35_)
        skele = skele.replace('_fuel_volume45_', _fuel_volume45_)

        skele = skele.replace('_matrix_density_', _matrix_density_)
        skele = skele.replace('_matrix_temperature_', _matrix_temperature_)

        skele = skele.replace('_core_density_', _core_density_)
        skele = skele.replace('_shell_density_', _shell_density_)
        skele = skele.replace('_coolant_density_', _coolant_density_)
        skele = skele.replace('_core_temperature_', _core_temperature_)
        skele = skele.replace('_shell_temperature_', _shell_temperature_)
        skele = skele.replace('_coolant_density_', _coolant_density_)

        skele = skele.replace(
            '_coolant_outlet_density_',
            _coolant_outlet_density_)
        skele = skele.replace('_coolant_density_', _coolant_density_)
        skele = skele.replace(
            '_coolant_inlet_density_',
            _coolant_inlet_density_)

        skele = skele.replace(
            '_inner_solid_reflector_density_',
            _inner_solid_reflector_density_)
        skele = skele.replace(
            '_density_inner_porous_reflector_',
            _density_inner_porous_reflector_)
        skele = skele.replace(
            '_density_outer_porous_reflector_',
            _density_outer_porous_reflector_)

        skele = skele.replace(
            '_coolant_outlet_temperature_',
            _coolant_outlet_temperature_)
        skele = skele.replace('_coolant_temperature_', _coolant_temperature_)
        skele = skele.replace(
            '_coolant_inlet_temperature_',
            _coolant_inlet_temperature_)

        skele = skele.replace('_kernel_radius_', _kernel_radius_)
        skele = skele.replace('_triso_hpitch_', _triso_hpitch_)

        skele = skele.replace('_pebble_hpitch1_', _pebble_hpitch1_)
        skele = skele.replace('_pebble_hpitch2_', _pebble_hpitch2_)
        skele = skele.replace('_pebble_core_radius_', _pebble_core_radius_)
        skele = skele.replace('_pebble_active_radius_', _pebble_active_radius_)
        skele = skele.replace('_pebble_shell_radius_', _pebble_shell_radius_)

        skele = skele.replace(
            '_inner_solid_reflector_radius_',
            _inner_solid_reflector_radius_)
        skele = skele.replace(
            '_outer_solid_refelctor_radius_',
            _outer_solid_refelctor_radius_)
        skele = skele.replace(
            '_inner_solid_reflector_interface_',
            _inner_solid_reflector_interface_)
        skele = skele.replace(
            '_outer_solid_refelctor_interface_',
            _outer_solid_refelctor_interface_)
        skele = skele.replace('_shield_radius_', _shield_radius_)
        skele = skele.replace('_core_barrel_radius_', _core_barrel_radius_)
        skele = skele.replace('_down_comer_radius_', _down_comer_radius_)
        skele = skele.replace(
            '_reactor_vessel_radius_',
            _reactor_vessel_radius_)

        skele = skele.replace(
            '_entrance_fuel_radius1_',
            _entrance_fuel_radius1_)
        skele = skele.replace(
            '_entrance_fuel_radius2_',
            _entrance_fuel_radius2_)
        skele = skele.replace(
            '_entrance_fuel_radius3_',
            _entrance_fuel_radius3_)
        skele = skele.replace(
            '_entrance_fuel_radius4_',
            _entrance_fuel_radius4_)
        skele = skele.replace(
            '_entrance_pebble_reflector_radius_',
            _entrance_pebble_reflector_radius_)

        skele = skele.replace('_expansion_fuel_cone1_', _expansion_fuel_cone1_)
        skele = skele.replace('_expansion_fuel_cone2_', _expansion_fuel_cone2_)
        skele = skele.replace('_expansion_fuel_cone3_', _expansion_fuel_cone3_)
        skele = skele.replace('_expansion_fuel_cone4_', _expansion_fuel_cone4_)
        skele = skele.replace(
            '_expansion_pebble_reflector_core_',
            _expansion_pebble_reflector_core_)

        skele = skele.replace('_active_fuel_radius1_', _active_fuel_radius1_)
        skele = skele.replace('_active_fuel_radius2_', _active_fuel_radius2_)
        skele = skele.replace('_active_fuel_radius3_', _active_fuel_radius3_)
        skele = skele.replace('_active_fuel_radius4_', _active_fuel_radius4_)
        skele = skele.replace(
            '_active_pebble_reflector_radius1_',
            _active_pebble_reflector_radius1_)
        skele = skele.replace(
            '_active_pebble_reflector_radius2_',
            _active_pebble_reflector_radius2_)

        skele = skele.replace(
            '_converging_fuel_cone1_',
            _converging_fuel_cone1_)
        skele = skele.replace(
            '_converging_fuel_cone2_',
            _converging_fuel_cone2_)
        skele = skele.replace(
            '_converging_fuel_cone3_',
            _converging_fuel_cone3_)
        skele = skele.replace(
            '_converging_fuel_cone4_',
            _converging_fuel_cone4_)
        skele = skele.replace(
            '_converging_pebble_reflector_cone_',
            _converging_pebble_reflector_cone_)

        skele = skele.replace(
            '_defueling_fuel_radius1_',
            _defueling_fuel_radius1_)
        skele = skele.replace(
            '_defueling_fuel_radius2_',
            _defueling_fuel_radius2_)
        skele = skele.replace(
            '_defueling_fuel_radius3_',
            _defueling_fuel_radius3_)
        skele = skele.replace(
            '_defueling_fuel_radius4_',
            _defueling_fuel_radius4_)
        skele = skele.replace(
            '_defueling_pebble_reflector_radius_',
            _defueling_pebble_reflector_radius_)

        skele = skele.replace(
            '_plenum_reactor_vessel_',
            _plenum_reactor_vessel_)
        skele = skele.replace(
            '_plenum_coolant_surface_',
            _plenum_coolant_surface_)
        skele = skele.replace('_entrance_bottom_', _entrance_bottom_)
        skele = skele.replace(
            '_entrance_axial_height_',
            _entrance_axial_height_)
        skele = skele.replace(
            '_expansion_axial_height_',
            _expansion_axial_height_)
        skele = skele.replace('_active_axial_height_', _active_axial_height_)
        skele = skele.replace('_bottom_control_rod_', _bottom_control_rod_)
        skele = skele.replace('_top_control_rod_', _top_control_rod_)
        skele = skele.replace(
            '_top_porous_inner_reflector_',
            _top_porous_inner_reflector_)
        skele = skele.replace(
            '_bottom_porous_outer_reflector_',
            _bottom_porous_outer_reflector_)
        skele = skele.replace(
            '_converging_axial_height_',
            _converging_axial_height_)
        skele = skele.replace(
            '_defueling_axial_height_',
            _defueling_axial_height_)
        skele = skele.replace(
            '_entrance2_axial_height_',
            _entrance2_axial_height_)
        skele = skele.replace('_repose_height_', _repose_height_)

        skele = skele.replace('_fuel_materials_', _fuel_materials_[:-1])
        skele = skele.replace(
            '_porous_graphite_inlet_',
            _porous_graphite_inlet_)
        skele = skele.replace(
            '_porous_graphite_outlet_',
            _porous_graphite_outlet_)
        skele = skele.replace('_graphite_inlet_', _graphite_inlet_)
        skele = skele.replace('_graphite_outlet_', _graphite_outlet_)
        skele = skele.replace('_graphite_pebble_core_', _graphite_pebble_core_)
        skele = skele.replace(
            '_graphite_pebble_shell_',
            _graphite_pebble_shell_)
        skele = skele.replace('_triso_matrix_', _triso_matrix_)
        skele = skele.replace(
            '_inlet_coolant_material_',
            _inlet_coolant_material_)
        skele = skele.replace(
            '_outlet_coolant_material_',
            _outlet_coolant_material_)
        skele = skele.replace('_coolant_material_', _coolant_material_)
        skele = skele.replace('_control_rod_density_', _control_rod_density_)
        skele = skele.replace(
            '_control_blade_density_',
            _control_blade_density_)
        skele = skele.replace('_control_rod_material_', _control_rod_material_)
        skele = skele.replace(
            '_control_blade_material_',
            _control_blade_material_)

        skele = skele.replace('_Heat_Exchangers_', _Heat_Exchangers_)
        skele = skele.replace(
            '_heat_exchanger_radius_',
            _heat_exchanger_radius_)
        skele = skele.replace(
            '_heat_exchanger_axial_height_',
            _heat_exchanger_axial_height_)
        skele = skele.replace(
            '_heat_exchanger_density_',
            _heat_exchanger_density_)

        skele = skele.replace('_tally_', _tally_)

        skele = skele.replace('_absorber_radius_', _absorber_radius_)
        skele = skele.replace('_absorber_clad_radius_', _absorber_clad_radius_)
        skele = skele.replace('_absorber_hpitch_', _absorber_hpitch_)
        skele = skele.replace('_absorber_pitch_', _absorber_pitch_)
        skele = skele.replace('_shield_inner_radius_', _shield_inner_radius_)
        skele = skele.replace('_absorber_density_', _absorber_density_)
        skele = skele.replace('_shield_clad_density_', _shield_clad_density_)
        skele = skele.replace('_shield_block_density_', _shield_block_density_)
        skele = skele.replace(
            '_shield_absorber_material_',
            _shield_absorber_material_)
        skele = skele.replace('_shield_clad_material_', _shield_clad_material_)

        skele = skele.replace('_wall_pf_', _wall_pf_)
        skele = skele.replace('_bulk_pf_', _bulk_pf_)

        self.inp = skele
        self.decf = skele.replace(
            '.70c', '').replace(
            '.71c', '').replace(
            '.72c', '').replace(
            '.73c', '').replace(
            'fill=1 ', ' ').replace(
            'fill=2 ', ' ').replace(
            '.70c', '').replace(
            '.71c', '').replace(
            '.72c', '').replace(
            '.73c', '').replace(
            'c CR withdrawn ', '').replace(
            'c CB withdrawn ', '')
        open('inp.1', 'w').write(skele)

    def BEAU(self, title=''):

        import os
        import math

        self.tallies = ['mocup']
        self.MCNP5()

        directory = 'Equilibrium_%s' % title

        rm = 'rm -r ../%s' % directory
        # print(rm)
        os.system(rm)

        copy = 'cp -r Equilibrium_BEAU2 ../%s' % directory
        # print(copy)
        os.system(copy)

        skele = open('skeleton_BEAU8').read()

        _power_ = '%1.5E' % self.power
        _u235_ = '%1.5E' % (self.fuel.TRISO.mat_kernel.comp['922350'])
        _u238_ = '%1.5E' % (self.fuel.TRISO.mat_kernel.comp['922380'])

        skele = skele.replace('_power_', _power_)
        skele = skele.replace('_u235_', _u235_)
        skele = skele.replace('_u238_', _u238_)

        BEAU_loc = '../%s/BEAU8.py' % directory
        open(BEAU_loc, 'w').write(skele)

        inp = self.inp.replace(
            'c CR withdrawn ',
            '').replace(
            'c CB withdrawn ',
            '')
        inp_loc = '../%s/inp1' % directory
        open(inp_loc, 'w').write(inp)

        run = open('skeleton_run').read()
        if title != '':
            run = run.replace('_title_', title)
        run_loc = '../%s/run' % directory
        open(run_loc, 'w').write(run)
