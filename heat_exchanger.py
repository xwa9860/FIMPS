#!/usr/bin/python


class Heat_Exchanger:

    def __init__(self):

        import math
        import mocup

        self.T_primary = 650
        self.mat_primary = mocup.material()
        self.mat_primary.comp['30060'] = 0.0002
        self.mat_primary.comp['30070'] = 1.9999
        self.mat_primary.comp['40090'] = 1.0000
        self.mat_primary.comp['90190'] = 4.0000
        self.mat_primary = self.mat_primary * \
            ((2.28 - .000488*(self.T_primary))/self.mat_primary.mass())

        self.T_secondary = 650
        self.mat_secondary = mocup.material()
        self.mat_secondary.comp['30060'] = 0.0174375
        self.mat_secondary.comp['30070'] = 0.2150625
        self.mat_secondary.comp['110230'] = 0.0575
        self.mat_secondary.comp['190390'] = 0.195846
        self.mat_secondary.comp['190400'] = 0.0000252
        self.mat_secondary.comp['190410'] = 0.014133
        self.mat_secondary.comp[' 90190'] = 0.5
        self.mat_secondary = self.mat_secondary * \
            ((2.530 - .00073*self.T_secondary)/self.mat_secondary.mass())

        self.mat_tube = mocup.material()
        self.mat_tube.comp['280580'] = 0.043837414
        self.mat_tube.comp['280600'] = 0.016323137
        self.mat_tube.comp['280610'] = 0.000697987
        self.mat_tube.comp['280620'] = 0.002189101
        self.mat_tube.comp['280640'] = 0.000540385
        self.mat_tube.comp['240500'] = 0.000324558
        self.mat_tube.comp['240520'] = 0.006018055
        self.mat_tube.comp['240530'] = 0.000669524
        self.mat_tube.comp['240540'] = 0.000163572
        self.mat_tube.comp['420920'] = 0.001549149
        self.mat_tube.comp['420940'] = 0.000945063
        self.mat_tube.comp['420950'] = 0.001609409
        self.mat_tube.comp['420960'] = 0.001668675
        self.mat_tube.comp['420970'] = 0.000945537
        self.mat_tube.comp['420980'] = 0.002364712
        self.mat_tube.comp['421000'] = 0.000924854
        self.mat_tube.comp['260540'] = 0.000286536
        self.mat_tube.comp['260560'] = 0.004369387
        self.mat_tube.comp['260570'] = 0.000102966
        self.mat_tube.comp['260580'] = 1.28788E-05
        self.mat_tube = self.mat_tube*(8.86/self.mat_tube.mass())
        # Hastalloy N

        self.dr_secondary_channel = 1.11
        self.dr_tube = 0.16
        self.PD = 1.2

    def volume(self):

        import math

        self.r_secondary_channel = self.dr_secondary_channel
        self.r_tube = self.r_secondary_channel + self.dr_tube
        self.hpitch = self.PD*self.r_tube

        self.da_channel = math.pi*self.r_secondary_channel**2.
        self.da_tube = math.pi*self.r_tube**2. - self.da_channel
        self.A = 12./2.*math.tan(math.pi/6.)*self.hpitch**2.
        self.da_primary = self.A - self.da_channel - self.da_tube

    def mat(self):
        self.volume()
        self.mat_homogenized = self.mat_primary * \
            (self.da_primary/self.A) + self.mat_tube*(self.da_tube/self.A) + self.mat_secondary*(self.da_channel/self.A)

        return self.mat_homogenized


