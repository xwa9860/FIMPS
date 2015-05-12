#!/usr/bin/python
from triso import TRISO


class Pebble:

    def __init__(self):

        import mocup

        self.TRISO = TRISO()
        self.dr_core = 0.577
        self.dr_shell = 0.1
        self.r_shell = 1.5
        self.dr_active = self.r_shell - self.dr_shell - self.dr_core

        self.mat_core = mocup.material()
        self.mat_core.comp = {'60120': (8.50688E-02/.6022)}
        self.mat_active = self.TRISO.homogenize()
        self.mat_shell = mocup.material()
        self.mat_shell.comp = {'60120': (8.72394E-02/.6022)}

        self.mat_core_min = mocup.material()
        self.mat_core_max = mocup.material()
        self.mat_core_min.comp = {'60120': (2.50688E-02/.6022)}
        self.mat_core_max.comp = {'60120': (8.72394E-02/.6022)}

        self.mat_coolant = mocup.material()
        self.mat_coolant.comp = {
            '30060': 3.99919E-07,
            '30070': 3.95956E-02,
            '40090': 1.97980E-02,
            '90190': 7.91919E-02}

        self.min_TRISO = 1.
        self.CHM_err = .01
        self.rho_err = .01
        self.max_PF = .4
        self.max_iterations = 20

        self.Pebble_PF = 0.60
        self.T_kernel = 0
        self.T_coolant = 650 + 273.15
        self.volume()

    def volume(self):

        import math

        self.dr_active = self.r_shell - self.dr_shell - self.dr_core
        self.r_core = self.dr_core
        self.r_active = self.r_shell - self.dr_shell

        self.dv_core = 4./3.*math.pi*math.pow(self.r_core, 3.)
        self.dv_active = 4./3.*math.pi * \
            (math.pow(self.r_active, 3.)-math.pow(self.r_core, 3.))
        self.dv_shell = 4./3.*math.pi * \
            (math.pow(self.r_shell, 3.) - math.pow(self.r_active, 3.))
        self.V = self.dv_core + self.dv_active + self.dv_shell
        self.dv_coolant = self.V/self.Pebble_PF*(1-self.Pebble_PF)

        # FCC half pitch
        self.h_pitch = math.pow(4*self.V/self.Pebble_PF, (1./3.))/2.
        # Simple Hexagonal PRISM unit cell pitch
        self.pitch = 2 * \
            math.pow(self.V/(self.Pebble_PF*12.*math.tan(math.pi/6.)), (1./3.))

    def CHM(self):

        self.mat_active = self.TRISO.homogenize()
        self.volume()

        mat = self.mat_core*self.dv_core + self.mat_active * \
            self.dv_active + self.mat_shell*self.dv_shell

        self.mat_pebble = (mat + self.mat_coolant*self.dv_coolant)*(1/self.V)

        return (mat.comp['60120']/mat.heavy_moles()), (mat.mass()/self.V)

    def min_dr_active(self):

        self.dr_active = 2*self.TRISO.r_matrix*self.min_TRISO
        if self.dr_active > (self.r_shell - self.dr_shell):
            self.dr_active = (self.r_shell - self.dr_shell)

        self.dr_core = self.r_shell - self.dr_active - self.dr_shell
        self.volume()

    def generate(self, CHM, rho):

        # This method generates a fuel design with minimum average delta T

        # determine feasilibility of using a fuel design with 40% TRISO packing
        # fraction

        self.dr_core = 0
        self.TRISO.setPF(self.max_PF)
        min_CHM, max_rho = self.CHM()

        self.mat_core.comp['60120'] = self.mat_core_max.moles()
        self.min_dr_active()

        max_CHM, rho1 = self.CHM()

        self.mat_core.comp['60120'] = self.mat_core_min.moles()
        CHM1, min_rho = self.CHM()

        if min_CHM < CHM and CHM < max_CHM and min_rho < rho and rho < max_rho:
            # use design sequence 1
            #print("design sequence 1")
            self.design_sequence1(CHM, rho)

        elif CHM < min_CHM or rho > max_rho:

            print('not feasible')

        else:
            # use design sequence 2
            #print("design sequence 2")
            self.design_sequence2(CHM, rho)

    def design_sequence1(self, CHM, rho):

        # This method assumes a maximum TRISO Packing Fraction and iteratively
        # perturbes dr_core and density of the core

        import math

        # solve for volume fraction of active region in pebble, f

        self.volume()
        self.TRISO.homogenize()

        g = self.dv_shell/self.V

        f = (self.mat_shell.comp['60120'] * g + self.mat_core.comp['60120']
             * (1 - g)) / (CHM * self.mat_active.heavy_moles() + self.mat_core.comp
                           ['60120'] - self.mat_active.comp['60120'])

        if (f + g) < 1 and f > 0:
            self.dr_core = math.pow((1. - f - g), (1./3.))*self.r_shell
            self.volume()

            # solve for density of core to impose density

            self.mat_core.comp['60120'] = (
                rho * self.V - self.mat_shell.mass() * self.dv_shell -
                self.mat_active.mass() * self.dv_active) / (
                self.dv_core * 12.)

            iteration = 1

        else:
            iteration = 1 + self.max_iterations

        CHM1, rho1 = self.CHM()

        while (math.fabs(CHM1 - CHM)/CHM > self.CHM_err or math.fabs(rho1 - rho)/rho > self.rho_err) and iteration < self.max_iterations:

            iteration += 1

            f = (self.mat_shell.comp['60120'] * g + self.mat_core.comp['60120']
                 * (1 - g)) / (CHM * self.mat_active.heavy_moles() +
                               self.mat_core.comp['60120'] - self.mat_active.comp['60120'])
            if (f + g) < 1:
                self.dr_core = math.pow((1. - f - g), (1./3.))*self.r_shell
                self.volume()

                self.mat_core.comp['60120'] = (
                    rho * self.V - self.mat_shell.mass() * self.dv_shell -
                    self.mat_active.mass() * self.dv_active) / (
                    self.dv_core * 12.)

            else:
                self.dr_core = 0
                self.volume()

            if self.mat_core.comp['60120'] < self.mat_core_min.comp['60120']:
                self.mat_core.comp['60120'] = self.mat_core_min.moles()

                self.design_sequence3(CHM, rho)

            elif self.mat_core.comp['60120'] > self.mat_core_max.comp['60120']:
                self.mat_core.comp['60120'] = self.mat_core_max.moles()

                self.deisgn_sequence3(CHM, rho)

            CHM1, rho1 = self.CHM()

        if math.fabs(CHM1 - CHM)/CHM < self.CHM_err and math.fabs(rho1 - rho)/rho < self.rho_err and self.dr_active >= (self.min_TRISO*self.TRISO.r_matrix*2):

            print('successfully generated fuel design with design sequence 1')

        else:
            self.design_sequence2(CHM, rho)

    def design_sequence2(self, CHM, rho):

        # This method assumes a minimum annular thickness geometry and perturbs
        # PF and density of the core to generate the fuel design

        #print('design sequence 2')

        import math

        self.CHM2PF(CHM)

        self.mat_core.comp['60120'] = (
            rho * self.V - self.mat_shell.mass() * self.dv_shell -
            self.mat_active.mass() * self.dv_active) / (
            self.dv_core * 12.)

        CHM1, rho1 = self.CHM()
        iteration = 1

        while (math.fabs(CHM1 - CHM)/CHM > self.CHM_err or math.fabs(rho1 - rho)/rho > self.rho_err) and iteration < self.max_iterations:

            self.CHM2PF(CHM)
            self.mat_core.comp['60120'] = (
                rho * self.V - self.mat_shell.mass() * self.dv_shell -
                self.mat_active.mass() * self.dv_active) / (
                self.dv_core * 12.)

            if self.mat_core.comp['60120'] > self.mat_core_max.comp['60120']:
                self.mat_core.comp['60120'] = self.mat_core_max.moles()
                self.design_sequence3(CHM, rho)

            elif self.mat_core.comp['60120'] < self.mat_core_min.comp['60120']:
                self.mat_core.comp['60120'] = self.mat_core_max.moles()
                self.design_sequence3(CHM, rho)

            CHM1, rho1 = self.CHM()
            iteration += 1

        if math.fabs(CHM1 - CHM)/CHM < self.CHM_err and math.fabs(rho1 - rho)/rho < self.rho_err:

            print('successfully generated fuel design with design sequence 2')

        else:

            print(
                'did not successfully generated fuel design with design sequence 2')

    def design_sequence3(self, CHM, rho):

        import math

        #print('design sequence3')

        CHM1, rho1 = self.CHM()

        g = self.dv_shell/self.V

        iteration = 1

        while (math.fabs(CHM1 - CHM)/CHM > self.CHM_err or math.fabs(rho1 - rho)/rho > self.rho_err) and iteration < self.max_iterations:

            self.CHM2PF2(CHM)

            f = (rho - g*self.mat_shell.mass() - (1-g)*self.mat_core.mass()) / \
                (self.mat_active.mass() - self.mat_core.mass())

            if (f + g) < 1 and f > 0:

                self.dr_core = math.pow((1. - f - g), (1./3.))*self.r_shell
                self.volume()

            else:
                #print('core radius')
                #print("fuck my life")

                self.min_dr_active()
                self.volume()
                self.CHM2PF2(CHM)

            if self.TRISO.PF > self.max_PF:

                print('packing fraction')
                print('fuck my life')

                print(self.TRISO.PF)

            CHM1, rho1 = self.CHM()

            iteration += 1

    def design_sequence_solid(self, CHM, rho):

        import math

        #print('solid pebble design sequence')

        self.dr_core = 1.e-30
        self.volume()

        CHM1, rho1 = self.CHM()

        iteration = 1

        while (math.fabs(CHM1 - CHM)/CHM > self.CHM_err or math.fabs(rho1 - rho)/rho > self.rho_err) and iteration < self.max_iterations:

            self.CHM2PF2(CHM)

            if self.TRISO.PF > self.max_PF:
                #print('packing fraction | fuck my life')
                # print(self.TRISO.PF)
                poop

            rho_triso = (
                self.TRISO.dv_kernel * self.TRISO.mat_kernel.mass() +
                self.TRISO.dv_buffer * self.TRISO.mat_buffer.mass() +
                self.TRISO.dv_iPyC * self.TRISO.mat_iPyC.mass() +
                self.TRISO.dv_SiC * self.TRISO.mat_SiC.mass() +
                self.TRISO.dv_oPyC * self.TRISO.mat_oPyC.mass()) / (
                self.TRISO.PF * self.TRISO.V)

            rho_matrix = (rho - (self.dv_shell/self.V)*self.mat_shell.mass() - (self.dv_active/self.V) *
                          (1 - (self.TRISO.dv_matrix/self.TRISO.V))*rho_triso)/((self.dv_active*self.TRISO.dv_matrix)/(self.V*self.TRISO.V))

            if rho_matrix > 1.75:
                rho_matrix = 1.75
            elif rho_matrix < 0.5:
                rho_matrix = 0.5

            self.TRISO.mat_matrix.comp = {}
            self.TRISO.mat_matrix.comp['60120'] = rho_matrix/12.

            CHM1, rho1 = self.CHM()

            iteration += 1

        if (math.fabs(CHM1 - CHM)/CHM > self.CHM_err or math.fabs(rho1 - rho)/rho > self.rho_err):
            poop

    def CHM2PF2(self, CHM):

        import math

        PF = {'min': 1.05*self.TRISO.PF, 'mid': 0, 'max': self.TRISO.PF}
        CHMs = {'mid': 0}

        self.TRISO.setPF(PF['min'])
        CHMs['min'], rho1 = self.CHM()

        self.TRISO.setPF(PF['max'])
        CHMs['max'], rho1 = self.CHM()

        min = 0
        max = 0

        iteration = 1

        while math.fabs(CHMs['mid']-CHM)/CHM > self.CHM_err and iteration < self.max_iterations:

            iteration += 1

            PF['mid'] = (PF['min'] - PF['max']) / \
                (CHMs['min'] - CHMs['max'])*(CHM - CHMs['max']) + PF['max']

            if PF['mid'] < PF['min']*.1:
                PF['mid'] = PF['min']*.1
            elif PF['mid'] > (self.max_PF - PF['min'])*.9 + PF['min']:
                PF['mid'] = (self.max_PF - PF['min'])*.9 + PF['min']
            elif (min > 6) or (max > 6):
                PF['mid'] = .5*(PF['max'] + PF['min'])

            self.TRISO.setPF(PF['mid'])
            CHMs['mid'], rho1 = self.CHM()

            if PF['mid'] > PF['max']:
                PF['min'] = PF['max']
                CHMs['min'] = CHMs['max']
                PF['max'] = PF['mid']
                CHMs['max'] = CHMs['mid']
            elif PF['mid'] < PF['min']:
                PF['max'] = PF['min']
                CHMs['max'] = PF['min']
                PF['min'] = PF['mid']
                CHMs['min'] = CHMs['mid']
            else:
                if CHM > CHMs['mid']:
                    PF['max'] = PF['mid']
                    CHMs['max'] = CHMs['mid']
                else:
                    PF['min'] = PF['mid']
                    CHMs['min'] = CHMs['mid']

    def CHM2PF(self, CHM):

        # This method determines the PF to impose a CHM holding everything else
        # in the fuel design constanti

        import math

        PF = {'min': .5*self.max_PF, 'mid': 0, 'max': self.max_PF}
        CHMs = {'mid': 0}

        self.TRISO.setPF(PF['min'])
        self.min_dr_active()
        CHMs['min'], rho1 = self.CHM()

        self.TRISO.setPF(PF['max'])
        self.min_dr_active()
        CHMs['max'], rho1 = self.CHM()

        min = 0
        max = 0

        iteration = 1

        while math.fabs(CHMs['mid']-CHM)/CHM > self.CHM_err and iteration < 15:

            iteration += 1

            PF['mid'] = (PF['min'] - PF['max']) / \
                (CHMs['min'] - CHMs['max'])*(CHM - CHMs['max']) + PF['max']

            # check reasonableness of PF['mid']

            if PF['mid'] < PF['min']*.1:
                PF['mid'] = PF['min']*.1
            elif PF['mid'] > (self.max_PF - PF['min'])*.9 + PF['min']:
                PF['mid'] = (self.max_PF - PF['min'])*.9 + PF['min']
            elif (min > 6) or (max > 6):
                PF['mid'] = .5*(PF['max'] + PF['min'])

            self.TRISO.setPF(PF['mid'])
            self.min_dr_active()
            CHMs['mid'], rho1 = self.CHM()

            if PF['mid'] > PF['max']:
                PF['min'] = PF['max']
                CHMs['min'] = CHMs['max']
                PF['max'] = PF['mid']
                CHMs['max'] = CHM['mid']
            elif PF['mid'] < PF['min']:
                PF['max'] = PF['min']
                CHMs['max'] = PF['min']
                PF['min'] = PF['mid']
                CHMs['min'] = CHMs['mid']
            else:
                if CHM > CHMs['mid']:
                    PF['max'] = PF['mid']
                    CHMs['max'] = CHMs['mid']
                elif CHM < CHMs['mid']:
                    PF['min'] = PF['mid']
                    CHMs['min'] = CHMs['mid']

    def conductivities(self):

        self.k_shell = 15.
        self.k_active = self.TRISO.conductivities()
        p = self.mat_core.mass()/self.mat_shell.mass()
        # core conductivity adjusted for porosity
        self.k_core = self.k_shell*(1-p)/(1+2*p)

    def HT(self, T_bulk, Re, Power_Density):

        import math

        # generate coolant thermal properties based on T_bulk (K)

        self.volume()
        self.conductivities()

        rho = 2282. - 0.49*(T_bulk - 273.15)
        k1 = 0.119*math.pow((469+273.15), 0.5) * \
            math.pow((33.04/(rho/1000)), .667)/math.pow((33.04/2.17), 1.167)
        k2 = 0.0005*T_bulk+32/33.04 - .34
        k_coolant = (k1 + k2)/2.

        viscosity = .116*.001*math.exp(3755/T_bulk)
        Pr = 2368*viscosity/k_coolant

        Nu = 2 + 1.1*math.pow(Re, .6)*math.pow(Pr, .33)
        h = Nu*k_coolant/(2*self.r_shell*.01)

        dt_coolant = 1/(h*4*math.pi*math.pow(self.r_shell*.01,
                                             2.))*Power_Density*1E+06*(4./3.*math.pi*math.pow(self.r_shell*.01,
                                                                                              3.))/self.Pebble_PF

        dt_shell = (1. / (self.r_active * .01) - 1. / (self.r_shell * .01)) / (4 *
                                                                               math.pi * self.k_shell) * Power_Density * 1E+06 * (4. / 3. * math.pi
                                                                                                                                  * math.pow(self.r_shell * .01, 3.)) / self.Pebble_PF
        at_shell = (1./(2.*self.k_shell)*(math.pow((self.r_shell*.01),
                                                   2.)-math.pow((self.r_active*.01),
                                                                2.))/(4./3.*math.pi*(math.pow((self.r_shell*.01),
                                                                                              3.)-math.pow((self.r_active*.01),
                                                                                                           3.))) - 1./(4*math.pi*self.r_shell*.01*self.k_shell))*Power_Density*1E+06*(4./3.*math.pi*math.pow(self.r_shell*.01,
                                                                                                                                                                                                             3.))/self.Pebble_PF

        q = (Power_Density * 1E+06
             * (4. / 3. * math.pi * math.pow(self.r_shell * .01, 3.)) /
             self.Pebble_PF) / (self.dv_active * 1e-6)
        dt_active = q/(6*self.k_active)*math.pow((self.r_active*.01),
                                                 2.)*(1.-math.pow((self.r_core/self.r_active),
                                                                  2.)) + q/(3.*self.k_active)*math.pow((self.r_core*.01),
                                                                                                       3.)/(self.r_active*.01)*(1.-self.r_active/self.r_core)

        at_active = q*math.pi/(45.*self.k_active)*(4.*math.pow((self.r_active*.01),
                                                               5) - 20.*math.pow((self.r_core*.01),
                                                                                 6.)/(self.r_active*.01) + 36.*math.pow((self.r_core*.01),
                                                                                                                        5.) - 20*math.pow((self.r_core*.01),
                                                                                                                                          3)*math.pow((self.r_active*.01),
                                                                                                                                                      2.))/(self.dv_active*1.e-06)

        q = (Power_Density * 1E+06
             * (4. / 3. * math.pi * math.pow(self.r_shell * .01, 3.)) /
             self.Pebble_PF) / (self.dv_active / self.TRISO.V)

        self.TRISO.HT(q)

        self.T_shell = T_bulk + dt_coolant + at_shell
        self.T_active = T_bulk + dt_coolant + dt_shell + at_active
        self.T_core = T_bulk + dt_coolant + dt_shell + dt_active
        self.T_kernel = T_bulk + dt_coolant + \
            dt_shell + at_active + self.TRISO.ave_T
        self.max_T = T_bulk + dt_coolant + \
            dt_shell + dt_active + self.TRISO.max_T

        return self.T_kernel, self.max_T

    def inputMCNP5(self, T_bulk, Re, power_density):

        import math
        import mocup

        # T_bulk        - bulk temperature of the coolant in Kelvin (K)
        # Re            - Reynolds number of coolant in PB-AHTR, baseline value ~ 1200
        # power_density - power density of the unit cell (MW/m3)
    # power of unit cell / volume of unit cell (including pebble core, shell
    # and coolant)
        self.TRISO.matrix()
        self.volume()

        self.HT(T_bulk, Re, power_density)
        with open('skeleton_core', 'r+') as skeletonFile:
            skeleton = skeletonFile.read()
        # start generating tokens to sub out in skeleton file

        __kernel_molar_density__ = '%1.5E' % (
            self.TRISO.mat_kernel.moles()*.60221415)
        __kernel_volume__ = '%1.5E' % (
            0.5*self.TRISO.dv_kernel/self.TRISO.V*self.dv_active)
        __kernel_temperature__ = '%1.5E' % (self.T_kernel*8.6173324e-11)

        __matrix_density__ = '%1.5E' % (self.TRISO.matrix().mass())
        __matrix_temperature__ = '%1.5E' % (self.T_active*8.6173324e-11)

        __pebble_core_molar_density__ = '%1.5E' % (
            self.mat_core.moles()*.60221415)
        __pebble_shell_molar_density__ = '%1.5E' % (
            self.mat_shell.moles()*.60221415)

        __kernel_radius__ = '%1.5E' % (self.TRISO.dr_kernel)
        __TRISO_hpitch__ = '%1.5E' % (math.pow((self.TRISO.V), (1./3.))/2.)

        __Pebble_hpitch__ = '%1.5E' % self.h_pitch
        __Pebble_core_radius__ = '%1.5E' % self.r_core
        __Pebble_active_radius__ = '%1.5E' % self.r_active
        __Pebble_shell_radius__ = '%1.5E' % self.r_shell

        __matrix_carbon_molar_density__ = '%1.5E' % (
            self.TRISO.mat_hmatrix.comp['60120']*.60221415)
        __matrix_silicon_molar_density__ = '%1.5E' % (
            self.TRISO.mat_hmatrix.comp['140280']*.60221415)

        XS = {
            '70c': [0, (.600-.5*300)],
            '71c': [((600-.5*300)), ((900-.5*300))],
            '72c': [((900-.5*300)), (.1200-(2500-1200)*.5)],
            '73c': [((.1200+-2500-1200)*.5), 5778.],
            }
        for xs in list(XS.keys()):

            if self.T_kernel > XS[xs][0] and self.T_kernel < XS[xs][1]:
                __XS_kernel__ = xs
            if self.max_T > XS[xs][0] and self.max_T < XS[xs][1]:
                __XS_pebble_core__ = xs
            if self.T_shell > XS[xs][0] and self.T_shell < XS[xs][1]:
                __XS_pebble_shell__ = xs
            if T_bulk > XS[xs][0] and T_bulk < XS[xs][1]:
                __XS_coolant__ = xs

        # replace tokens with specific values

        skeleton = skeleton.replace(
            '__kernel_molar_density__',
            __kernel_molar_density__)
        skeleton = skeleton.replace('__kernel_volume__', __kernel_volume__)
        skeleton = skeleton.replace(
            '__kernel_temperature__',
            __kernel_temperature__)
        skeleton = skeleton.replace('__matrix_density__', __matrix_density__)
        skeleton = skeleton.replace(
            '__matrix_temperature__',
            __matrix_temperature__)
        skeleton = skeleton.replace(
            '__pebble_core_molar_density__',
            __pebble_core_molar_density__)
        skeleton = skeleton.replace(
            '__pebble_shell_molar_density__',
            __pebble_shell_molar_density__)
        skeleton = skeleton.replace('__kernel_radius__', __kernel_radius__)
        skeleton = skeleton.replace('__TRISO_hpitch__', __TRISO_hpitch__)
        skeleton = skeleton.replace('__Pebble_hpitch__', __Pebble_hpitch__)
        skeleton = skeleton.replace(
            '__Pebble_core_radius__',
            __Pebble_core_radius__)
        skeleton = skeleton.replace(
            '__Pebble_active_radius__',
            __Pebble_active_radius__)
        skeleton = skeleton.replace(
            '__Pebble_shell_radius__',
            __Pebble_shell_radius__)
        skeleton = skeleton.replace('__XS_kernel__', __XS_kernel__)
        skeleton = skeleton.replace('__XS_pebble_core__', __XS_pebble_core__)
        skeleton = skeleton.replace('__XS_pebble_shell__', __XS_pebble_shell__)
        skeleton = skeleton.replace('__XS_coolant__', __XS_coolant__)
        skeleton = skeleton.replace(
            '__matrix_carbon_molar_density__',
            __matrix_carbon_molar_density__)
        skeleton = skeleton.replace(
            '__matrix_silicon_molar_density__',
            __matrix_silicon_molar_density__)
        self.inp = skeleton

        self.inp_source = skeleton.replace(
            'kcode  10000  1.0 10 110',
            'ksrc 0 0 0\nkcode 1000 1 10 110')

        skeleton1 = 'unit cell depletion of PB-AHTR pebbles\nc ------------------------------------------------------------------------------\nc Cells\nc ------------------------------------------------------------------------------\nc\nc TRISO unit cell\nc\n1 1  __kernel_molar_density__ -1    tmp=__kernel_temperature__ imp:n=1 u=100 $ Fuel Kernel\n          vol=__kernel_volume__\n4 3 -__matrix_density__  1    tmp=__matrix_temperature__ imp:n=1 u=100 $ Silicon Carbide Layer\nc\n10 0       10 -11 12 -13 14 -15 lat=1 imp:n=1 fill=100 u=10 $ infinite lattice of TRISO particles\nc\nc Pebble unit cell\nc\n100 2 __pebble_core_molar_density__  -100     tmp=9.38858e-08        imp:n=1      $ Low Density Graphite Pebble Core\n101 0              -101 100                fill=10 imp:n=1      $ Active Region of Pebble\n102 4  __pebble_shell_molar_density__ -102 101 tmp=8.27695-08         imp:n=1      $ Pebble Shell\n103 5  -__coolant_density__ 102 110 -111 112 -113 114 -115 116 -117     $ Flibe Coolant\n                     tmp=7.95552E-08         imp:n=1 \nc\nc Boundary Conditions\nc\n999 0               -110:111:-112:113:-114:115:-116:117 imp:n=1 $ Vacuum Boundary Conditions\nc neutrons should be reflected and these Vacuum Boundary Conditions will not be implemented\n\nc ------------------------------------------------------------------------------\nc Surfaces\nc ------------------------------------------------------------------------------\nc\nc Surfaces of TRISO layers\nc\n1 so __kernel_radius__  $ Outer Radius of Kernel\nc\nc Boundaries of TRISO unit cell\nc\n10 px -__TRISO_hpitch__\n11 px  __TRISO_hpitch__\n12 py -__TRISO_hpitch__\n13 py  __TRISO_hpitch__\n14 pz -__TRISO_hpitch__\n15 pz  __TRISO_hpitch__\nc\nc Radial Surfaces of Pebble \nc\n100 so __Pebble_core_radius__ $ Outer Radius of Low Density Graphite Pebble Core\n101 so __Pebble_active_radius__       $ Outer Radius of Pebble Active Region\n102 so __Pebble_shell_radius__       $ Outer Radius of Pebble Shell\nc \nc Boundaries of Pebble unit cell\nc\n*110 px                  -__Pebble_hpitch__\n*111 px                   __Pebble_hpitch__\n*112 p   1 1.732050808 0 -__Pebble_pitch__\n*113 p   1 1.732050808 0  __Pebble_pitch__\n*114 p  -1 1.732050808 0 -__Pebble_pitch__\n*115 p  -1 1.732050808 0  __Pebble_pitch__\n*116 pz                  -__Pebble_hpitch__\n*117 pz                   __Pebble_hpitch__\n\nc ------------------------------------------------------------------------------\nc Data\nc ------------------------------------------------------------------------------\n__kernel_fuel__\nm2         6000.72c 1 $ Graphite (1200K)\nmt2       grph.65t\n__matrix_mat__\nm4         6000.72c 1 $ Pebble Shell\nmt4       grph.64t\n__coolant_mat__\n__single_mat__\n__tallies__\nkcode 10000 1 10 110 \nksrc 0 0 0\n\n'

        mat = mocup.material()
        mat.addnux()
        self.TRISO.mat_kernel = self.TRISO.mat_kernel + mat

        __kernel_volume__ = '%1.5E' % (
            self.dv_active*self.TRISO.dv_kernel/self.TRISO.V)
        __coolant_density__ = '%1.5E' % self.mat_coolant.mass()
        __Pebble_hpitch__ = '%1.5E' % (self.pitch/2.)
        __Pebble_pitch__ = '%1.5E' % self.pitch

        self.TRISO.mat_kernel.mocup_strings(lib=__XS_kernel__)

        skeleton1 = skeleton1.replace(
            '__kernel_molar_density__',
            __kernel_molar_density__)
        skeleton1 = skeleton1.replace(
            '__kernel_temperature__',
            __kernel_temperature__)
        skeleton1 = skeleton1.replace('__kernel_volume__', __kernel_volume__)
        skeleton1 = skeleton1.replace('__matrix_density__', __matrix_density__)
        skeleton1 = skeleton1.replace(
            '__matrix_temperature__',
            __matrix_temperature__)
        skeleton1 = skeleton1.replace(
            '__pebble_core_molar_density__',
            __pebble_core_molar_density__)
        skeleton1 = skeleton1.replace(
            '__pebble_shell_molar_density__',
            __pebble_shell_molar_density__)
        skeleton1 = skeleton1.replace(
            '__coolant_density__',
            __coolant_density__)
        skeleton1 = skeleton1.replace('__kernel_radius__', __kernel_radius__)
        skeleton1 = skeleton1.replace('__TRISO_hpitch__', __TRISO_hpitch__)

        skeleton1 = skeleton1.replace(
            '__Pebble_core_radius__',
            __Pebble_core_radius__)
        skeleton1 = skeleton1.replace(
            '__Pebble_active_radius__',
            __Pebble_active_radius__)
        skeleton1 = skeleton1.replace(
            '__Pebble_shell_radius__',
            __Pebble_shell_radius__)
        skeleton1 = skeleton1.replace('__Pebble_hpitch__', __Pebble_hpitch__)
        skeleton1 = skeleton1.replace('__Pebble_pitch__', __Pebble_pitch__)
        skeleton1 = skeleton1.replace(
            '__kernel_fuel__',
            self.TRISO.mat_kernel.mcf(
                1,
                __XS_kernel__))
        skeleton1 = skeleton1.replace(
            '__matrix_mat__',
            self.TRISO.mat_hmatrix.mcf(
                3,
                __XS_kernel__))
        skeleton1 = skeleton1.replace(
            '__coolant_mat__',
            self.mat_coolant.mcf(
                5,
                __XS_coolant__))
        # print(dir(self.TRISO.mat_kernel))
        skeleton1 = skeleton1.replace(
            '__single_mat__',
            self.TRISO.mat_kernel.single_mat)
        skeleton1 = skeleton1.replace(
            '__tallies__',
            self.TRISO.mat_kernel.tally)

        skeleton1 = skeleton1.replace('_cells_', '1')

        self.SHC_inp = skeleton1

    def mocup(
            self,
            title='default',
            max_BU=5.e5,
            dBU1=200.,
            dBU2=1.e4,
            nBU1=5,
            T_bulk=650 +
            273.15,
            Re=1200,
            power_density=20):

        # This method sets up a mocup calculation

        import time
        import os
        import math

        if title == 'default':
            title = time.strftime("%a_%d_%b_%Y_%H:%M")

        self.inputMCNP5(T_bulk, Re, power_density)

        BU_vector = []

        for i in range(nBU1):
            BU_vector.append(dBU1)

        max_BU = max_BU - nBU1*dBU1
        nBU2 = math.ceil(max_BU/dBU2)

        for i in range(nBU2):
            BU_vector.append(dBU2)

        _power_vector_ = ''
        _time_vector_ = ''

        for dBU in BU_vector:
            _power_vector_ += '%1.5E,' % (power_density *
                                          self.V/self.Pebble_PF*1e-6)
            _time_vector_ += '%1.5E,' % (dBU*self.TRISO.mat_kernel.mass()*self.TRISO.dv_kernel *
                                         self.dv_active/(power_density*self.TRISO.V*self.V/self.Pebble_PF))

        skeleton = open('skeleton').read()
        run_mocup = "#! /bin/sh\n#\n#$ -N mocup\n#$ -cwd\n#$ -pe ompi 8\n#$ -S /bin/bash\n#$ -q x.q\n#$ -V\n\nmpiexec -x LD_LIBRARY_PATH  mcnp5.mpi i=inp.1 o=outp.1 mc=mctal.1 runtpe=runtp1\nmv srctp source\nrm outp.1 mctal.1 runtp1\nsed s/'ksrc 0 0 0'/'c'/ inp.1 >> tmp\nmv tmp inp.1\npython gofiss.py\n"

        skeleton = skeleton.replace('_title_', title)
        skeleton = skeleton.replace('_power_vector_', _power_vector_)
        skeleton = skeleton.replace('_time_vector_', _time_vector_)

        dir = '../depletion_%s' % title
        inp_loc = dir + '/inp.1'
        mocup_loc = '%s/mocup.py' % dir
        gofiss_loc = '%s/gofiss.py' % dir
        run_mocup_loc = '%s/run_mocup.sh' % dir
        moi_file = '%s/moi_files' % dir

        mkdir = '\\rm %s -r' % dir
        # print(mkdir)
        os.system(mkdir)

        mkdir = 'mkdir %s' % dir
        # print(mkdir)
        os.system(mkdir)

        mkdir = 'mkdir %s' % moi_file
        # print(mkdir)
        os.system(mkdir)

        open(inp_loc, 'w').write(self.SHC_inp)
        open(run_mocup_loc, 'w').write(run_mocup)
        open(gofiss_loc, 'w').write(skeleton)

        copy = '\cp Equilibrium_BEAU/mocup.py %s' % (mocup_loc)
        # print(copy)
        os.system(copy)

    def BEAU(self, title, T_bulk, Re, power_density):

        import os

        dir = '~/Equilibrium_%s' % title

        inp_loc = '\cp inp.1 ' + dir + '/inp.1'
        BEAU_loc = '\cp BEAU7.py ' + dir + '/BEAU7.py'
        source_loc = '\cp inp ' + dir + '/inp'

        copy = '\\rm -r %s' % dir
        # print(copy)
        os.system(copy)

        copy = '\cp -r Equilibrium_BEAU %s' % dir
        # print(copy)
        os.system(copy)

        self.inputMCNP5(T_bulk, Re, power_density)

        skeleton = "#!/usr/bin/env python\n\nimport BEAU\nimport mocup\n\n# -----------------------------------------------------------------------------\n# DEFINE EQUILIBRIUM DEPLETION ANALYSIS PROBLEM\n# -----------------------------------------------------------------------------\n\n# initialize Equilibrium Depletion\nEq_Cycle = BEAU.equilibrium_cycle()\n\n# DEFINE POWER OF SYSTEM IN MWth\nEq_Cycle.power = (__unit_cell_power__)\n        # 8 1/2 pebbles with a volume of 1.17810E+01 cc with a power desnity of 30 MW/m3\n\n# DEFINE LOADING PATTER\n        # 'Y':'X' EOEC material from X is advanced to cell Y\nEq_Cycle.loading_pattern = {\n                                                        '20':'10',\n                                                        '30':'20',\n                                                        '40':'30',\n                                                        '50':'40',\n                                                        '60':'50',\n                                                        '70':'60',\n                                                        '80':'70',\n                                                        }\n\n# DEFINE THE BURNUPS OF ALL FUEL PROGRESSIONS\nEq_Cycle.burnup = {'seed': 1.5E+5}\n\n# DEFINE THE MAKE UP MATERIAL FOR EACH PROGRESSION\n\n# initiate material object\nmat = mocup.material()\n# define its composition vector\nmat.comp['922350'] = __seed_u235_molar_mass__\nmat.comp['922380'] = __seed_u238_molar_mass__\n\n# link material object to the progression key\nEq_Cycle.makeup = {'seed':mat}\n\n# INITIATE A PROGRESSION\nprogression = Eq_Cycle.fuel_progression()\n\n# DEFINE THE ORDER IN WHICH FUEL IS ADVANCED\nprogression.cells = ['10','20','30','40','50','60','70','80']\n\n# DEFINE THE MAKEUP FUEL VECTOR FOR THE PROGRESSION\nprogression.feed = Eq_Cycle.makeup['seed']\n\n# DEFINE THE STANDARD ORIGEN LIBRARY TO BE USED\nprogression.library = 'pwru50'\n\n# LINK PROGRESSION OBJECTS TO KEYS\nEq_Cycle.progressions = {'seed':progression,}\n\n# -----------------------------------------------------------------------------\n# INITIATE EQUILIBRIUM DEPLETION ANALYSIS\n# -----------------------------------------------------------------------------\nEq_Cycle.dBU1 = 1000\nEq_Cycle.dBU2 = 10000\nEq_Cycle.nBU1 = 3\n\nEq_Cycle.search_keff(2.2E+05,2.8e+05, 1.000)"
        BOC_fuel = self.TRISO.mat_kernel * \
            (self.TRISO.dv_kernel/self.TRISO.V*self.dv_active*(1./2.))

        __unit_cell_power__ = '%1.5E' % (
            self.V/self.Pebble_PF*4.*power_density*1e-6)
        __seed_u235_molar_mass__ = '%1.5E' % BOC_fuel.comp['922350']
        __seed_u238_molar_mass__ = '%1.5E' % BOC_fuel.comp['922380']

        skeleton = skeleton.replace('__unit_cell_power__', __unit_cell_power__)
        skeleton = skeleton.replace(
            '__seed_u235_molar_mass__',
            __seed_u235_molar_mass__)
        skeleton = skeleton.replace(
            '__seed_u238_molar_mass__',
            __seed_u238_molar_mass__)

        open('inp.1', 'w').write(self.inp)
        open('BEAU7.py', 'w').write(skeleton)
        open('inp', 'w').write(self.inp_source)

        # print(inp_loc)
        os.system(inp_loc)
        # print(BEAU_loc)
        os.system(BEAU_loc)
        # print(source_loc)
        os.system(source_loc)

    def SERPENT(self):

        self.CHM()

        # This method makes a single pebble unit cell SERPENT input deck

        skeleton = '% --- FHR TRISO Unit Cell ------------------------------------------------------\n\nset title "FHR TRISO Unit Cell"\n\n% --- surfaces -----------------------------------------------------------------\n\nsurf 1 sph  0 0 0 _r_kernel_       % fuel kernel\nsurf 2 sph  0 0 0 _r_buffer_       % buffer\nsurf 3 sph  0 0 0 _r_iPyC_         % iPyC\nsurf 4 sph  0 0 0 _r_SiC_          % SiC\nsurf 5 sph  0 0 0 _r_oPyC_         % oPyC\nsurf 6 cube 0 0 0 _r_TRISO_hpitch_ % unit cell boundaries\n\nsurf 10 sph       0 0 0 _r_core_\nsurf 11 sph       0 0 0 _r_active_\nsurf 12 sph       0 0 0 _r_shell_\nsurf 13 hexxprism 0 0 _pebble_hexpitch_ -_pebble_hexpitch_ _pebble_hexpitch_ \n\n% --- cells --------------------------------------------------------------------\n\nlat 1 6 0 0 _r_TRISO_pitch_ 2\n\n% triso unit cell\n\ncell 1 2 fuel   -1   % fuel kernel\ncell 2 2 buffer -2 1 % buffer\ncell 3 2 PyC    -3 2 % iPyC\ncell 4 2 SiC    -4 3 % SiC\ncell 5 2 PyC    -5 4 % oPyC\ncell 6 2 matrix    5 % matrix\n\n% pebble unit cell\n\ncell 10 0 low_density_graphite -10    % low density kernel\ncell 11 0 fill 1               -11 10 % active region\ncell 12 0 shell                -12 11 % pebble shell\ncell 13 0 coolant              -13 12 % coolant\ncell 99 0 outside                  13 % reflective boundary conditions\n\n% --- boundary conditions ------------------------------------------------------\n\nset bc 3 % reflective boundary conditions\n\n% --- materials ----------------------------------------------------------------\n\n_mat_fuel_\n\n_mat_buffer_\n\n_mat_PyC_\n\n_mat_SiC_\n\n_mat_matrix_\n\n_mat_pebble_core_\n\n_mat_shell_\n\n_mat_coolant_\n\nset ures 1 2 1e-09\n 92235.06c 92238.06c\n% --- thermal scattering data --------------------------------------------------\n\ntherm buffer               grj2.18t\ntherm PyC                  grj2.18t\ntherm matrix               grj2.18t\ntherm low_density_graphite grj2.18t\ntherm shell                grj2.18t\n\n% --- Cross Section Library file path ------------------------------------------\n\nset acelib "/usr/local/SERPENT/xsdata/endfb7/sss_endfb7u.xsdata"\n\n% --- tallies ------------------------------------------------------------------\n\ndet "kernel spectrum" de 238g dm fuel dv 6.54498E-05\nene 238g 1\n     1.0000E-11 1.0000E-10 5.0000E-10 7.5000E-10 1.0000E-09 1.2000E-09\n     1.5000E-09 2.0000E-09 2.5000E-09 3.0000E-09 4.0000E-09 5.0000E-09\n     7.5000E-09 1.0000E-08 2.5300E-08 3.0000E-08 4.0000E-08 5.0000E-08\n     6.0000E-08 7.0000E-08 8.0000E-08 9.0000E-08 1.0000E-07 1.2500E-07\n     1.5000E-07 1.7500E-07 2.0000E-07 2.2500E-07 2.5000E-07 2.7500E-07\n     3.0000E-07 3.2500E-07 3.5000E-07 3.7500E-07 4.0000E-07 4.5000E-07\n     5.0000E-07 5.5000E-07 6.0000E-07 6.2500E-07 6.5000E-07 7.0000E-07\n     7.5000E-07 8.0000E-07 8.5000E-07 9.0000E-07 9.2500E-07 9.5000E-07\n     9.7500E-07 1.0000E-06 1.0100E-06 1.0200E-06 1.0300E-06 1.0400E-06\n     1.0500E-06 1.0600E-06 1.0700E-06 1.0800E-06 1.0900E-06 1.1000E-06\n     1.1100E-06 1.1200E-06 1.1300E-06 1.1400E-06 1.1500E-06 1.1750E-06\n     1.2000E-06 1.2250E-06 1.2500E-06 1.3000E-06 1.3500E-06 1.4000E-06\n     1.4500E-06 1.5000E-06 1.5900E-06 1.6800E-06 1.7700E-06 1.8600E-06\n     1.9400E-06 2.0000E-06 2.1200E-06 2.2100E-06 2.3000E-06 2.3800E-06\n     2.4700E-06 2.5700E-06 2.6700E-06 2.7700E-06 2.8700E-06 2.9700E-06\n     3.0000E-06 3.0500E-06 3.1500E-06 3.5000E-06 3.7300E-06 4.0000E-06\n     4.7500E-06 5.0000E-06 5.4000E-06 6.0000E-06 6.2500E-06 6.5000E-06\n     6.7500E-06 7.0000E-06 7.1500E-06 8.1000E-06 9.1000E-06 1.0000E-05\n     1.1500E-05 1.1900E-05 1.2900E-05 1.3750E-05 1.4400E-05 1.5100E-05\n     1.6000E-05 1.7000E-05 1.8500E-05 1.9000E-05 2.0000E-05 2.1000E-05\n     2.2500E-05 2.5000E-05 2.7500E-05 3.0000E-05 3.1250E-05 3.1750E-05\n     3.3250E-05 3.3750E-05 3.4600E-05 3.5500E-05 3.7000E-05 3.8000E-05\n     3.9100E-05 3.9600E-05 4.1000E-05 4.2400E-05 4.4000E-05 4.5200E-05\n     4.7000E-05 4.8300E-05 4.9200E-05 5.0600E-05 5.2000E-05 5.3400E-05\n     5.9000E-05 6.1000E-05 6.5000E-05 6.7500E-05 7.2000E-05 7.6000E-05\n     8.0000E-05 8.2000E-05 9.0000E-05 1.0000E-04 1.0800E-04 1.1500E-04\n     1.1900E-04 1.2200E-04 1.8600E-04 1.9250E-04 2.0750E-04 2.1000E-04\n     2.4000E-04 2.8500E-04 3.0500E-04 5.5000E-04 6.7000E-04 6.8300E-04\n     9.5000E-04 1.1500E-03 1.5000E-03 1.5500E-03 1.8000E-03 2.2000E-03\n     2.2900E-03 2.5800E-03 3.0000E-03 3.7400E-03 3.9000E-03 6.0000E-03\n     8.0300E-03 9.5000E-03 1.3000E-02 1.7000E-02 2.5000E-02 3.0000E-02\n     4.5000E-02 5.0000E-02 5.2000E-02 6.0000E-02 7.3000E-02 7.5000E-02\n     8.2000E-02 8.5000E-02 1.0000E-01 1.2830E-01 1.5000E-01 2.0000E-01\n     2.7000E-01 3.3000E-01 4.0000E-01 4.2000E-01 4.4000E-01 4.7000E-01\n     4.9952E-01 5.5000E-01 5.7300E-01 6.0000E-01 6.7000E-01 6.7900E-01\n     7.5000E-01 8.2000E-01 8.6110E-01 8.7500E-01 9.0000E-01 9.2000E-01\n     1.0100E+00 1.1000E+00 1.2000E+00 1.2500E+00 1.3170E+00 1.3560E+00\n     1.4000E+00 1.5000E+00 1.8500E+00 2.3540E+00 2.4790E+00 3.0000E+00\n     4.3040E+00 4.8000E+00 6.4340E+00 8.1873E+00 1.0000E+01 1.2840E+01\n     1.3840E+01 1.4550E+01 1.5683E+01 1.7333E+01 2.0000E+01\n\n% --- Neutron Population data --------------------------------------------------\n\nset pop 10000 100 100 1.0\n\nplot 3 1000 1000 \n\n\n\n\n'

        _r_kernel_ = '%1.5E' % self.TRISO.r_kernel
        _r_buffer_ = '%1.5E' % self.TRISO.r_buffer
        _r_iPyC_ = '%1.5E' % self.TRISO.r_iPyC
        _r_SiC_ = '%1.5E' % self.TRISO.r_SiC
        _r_oPyC_ = '%1.5E' % self.TRISO.r_oPyC
        _r_TRISO_hpitch_ = '%1.5E' % self.TRISO.hpitch
        _r_TRISO_pitch_ = '%1.5E' % (2*self.TRISO.hpitch)

        _r_core_ = '%1.5E' % self.r_core
        _r_active_ = '%1.5E' % self.r_active
        _r_shell_ = '%1.5E' % self.r_shell
        _pebble_hexpitch_ = '%1.5E' % (self.pitch/2.)

        if self.T_kernel == 0:
            print(
                'you must run the HT method of the pebble object to get the temperature distribution in the pebble')

        self.TRISO.mat_kernel.scf(
            'fuel',
            temperature=self.T_kernel,
            library='06c')
        self.TRISO.mat_buffer.scf(
            'buffer',
            temperature=self.T_active,
            library='06c')
        self.TRISO.mat_iPyC.scf('PyC', temperature=self.T_active, library='06c')
        self.TRISO.mat_SiC.scf('SiC', temperature=self.T_active, library='06c')
        self.TRISO.mat_matrix.scf(
            'matrix',
            temperature=self.T_active,
            library='06c')

        self.mat_core.scf(
            'low_density_graphite',
            temperature=self.T_core,
            library='06c')
        self.mat_shell.scf('shell', temperature=self.T_shell, library='06c')
        self.mat_coolant.scf(
            'coolant',
            temperature=self.T_coolant,
            library='06c')

        # print(self.TRISO.mat_kernel.mat_card)

        _mat_fuel_ = self.TRISO.mat_kernel.mat_card
        _mat_buffer_ = self.TRISO.mat_buffer.mat_card
        _mat_PyC_ = self.TRISO.mat_iPyC.mat_card
        _mat_SiC_ = self.TRISO.mat_SiC.mat_card
        _mat_matrix_ = self.TRISO.mat_matrix.mat_card
        _mat_pebble_core_ = self.mat_core.mat_card
        _mat_shell_ = self.mat_shell.mat_card
        _mat_coolant_ = self.mat_coolant.mat_card

        skeleton = skeleton.replace('_r_buffer_', _r_buffer_)
        skeleton = skeleton.replace('_r_iPyC_', _r_iPyC_)
        skeleton = skeleton.replace('_r_SiC_', _r_SiC_)
        skeleton = skeleton.replace('_r_oPyC_', _r_oPyC_)
        skeleton = skeleton.replace('_r_TRISO_hpitch_', _r_TRISO_hpitch_)

        skeleton = skeleton.replace('_r_kernel_', _r_kernel_)
        skeleton = skeleton.replace('_r_TRISO_pitch_', _r_TRISO_pitch_)
        skeleton = skeleton.replace('_r_core_', _r_core_)
        skeleton = skeleton.replace('_r_active_', _r_active_)
        skeleton = skeleton.replace('_r_shell_', _r_shell_)
        skeleton = skeleton.replace('_pebble_hexpitch_', _pebble_hexpitch_)

        skeleton = skeleton.replace('_mat_fuel_', _mat_fuel_)
        skeleton = skeleton.replace('_mat_buffer_', _mat_buffer_)
        skeleton = skeleton.replace('_mat_PyC_', _mat_PyC_)
        skeleton = skeleton.replace('_mat_SiC_', _mat_SiC_)
        skeleton = skeleton.replace('_mat_matrix_', _mat_matrix_)

        skeleton = skeleton.replace('_mat_pebble_core_', _mat_pebble_core_)
        skeleton = skeleton.replace('_mat_shell_', _mat_shell_)
        skeleton = skeleton.replace('_mat_coolant_', _mat_coolant_)

        self.serpent = skeleton


