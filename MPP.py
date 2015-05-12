import FIMPS

# for power_density in [23.3,28.3,33.3]:
for power_density in [23]:
    for active_thickness in [70]:
        # for pebble_reflector_thickness in [13,20,30]:
        for pebble_reflector_thickness in [20]:
            if (active_thickness + pebble_reflector_thickness + 10) > (70 + 20 + 40):
                self.r_reactor_vessel += (active_thickness +
                                          pebble_reflector_thickness +
                                          10) - (70 + 20 + 40)
            #case_ID = 'MPP_PD_%d_AT_%d_PRT_%d' % (power_density,active_thickness, pebble_reflector_thickness)
            case_ID = 'Mk1_Simple'
            print('case ID: %s' % case_ID)

            dr_total = 70 + 20 + 40

            PBFHR = FIMPS.PBFHR45()
            PBFHR.power = 236.
            PBFHR.dr_active_fuel = active_thickness
            PBFHR.dr_active_pebble_reflector = pebble_reflector_thickness
            PBFHR.dr_active_solid_reflector = dr_total - \
                active_thickness - pebble_reflector_thickness

            PBFHR.passes = 8
            PBFHR.control_rods = 8
            PBFHR.control_blades = 8
            PBFHR.r_control_rods = 6.5
            PBFHR.wall_PF = .57
            PBFHR.bulk_PF = .61
            PBFHR.fuel_power_density = power_density
            PBFHR.tallies = ['economy', 'spectrum']
            PBFHR.generate()
            # PBFHR.BEAU(title=case_ID)
            PBFHR.MCNP5()
#			print(PBFHR.decf)
