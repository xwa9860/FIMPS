#!/usr/bin/python


class shield:

    def __init__(
            self,
            absorber_density=2.52,
            clad_density=8.0,
            matrix_density=1.74):

        import mocup

        # define material composition of absorber
        self.mat_absorber = mocup.material()
        self.mat_absorber.comp = {
            '50100': 3.64662E-02,
            '50110': 1.46139e-01,
            '60120': 4.56512E-02}
        self.mat_absorber = self.mat_absorber * \
            (absorber_density/self.mat_absorber.mass())
        self.absorber_density = self.mat_absorber.mass()

        self.mat_matrix = mocup.material()
        self.mat_matrix.comp = {'60120': 1.45e-01}
        self.mat_matrix = self.mat_matrix * \
            (matrix_density/self.mat_matrix.mass())

        self.mat_clad = mocup.material()
        self.mat_clad.comp['280580'] = 0.043837414
        self.mat_clad.comp['280600'] = 0.016323137
        self.mat_clad.comp['280610'] = 0.000697987
        self.mat_clad.comp['280620'] = 0.002189101
        self.mat_clad.comp['280640'] = 0.000540385
        self.mat_clad.comp['240500'] = 0.000324558
        self.mat_clad.comp['240520'] = 0.006018055
        self.mat_clad.comp['240530'] = 0.000669524
        self.mat_clad.comp['240540'] = 0.000163572
        self.mat_clad.comp['420920'] = 0.001549149
        self.mat_clad.comp['420940'] = 0.000945063
        self.mat_clad.comp['420950'] = 0.001609409
        self.mat_clad.comp['420960'] = 0.001668675
        self.mat_clad.comp['420970'] = 0.000945537
        self.mat_clad.comp['420980'] = 0.002364712
        self.mat_clad.comp['421000'] = 0.000924854
        self.mat_clad.comp['260540'] = 0.000286536
        self.mat_clad.comp['260560'] = 0.004369387
        self.mat_clad.comp['260570'] = 0.000102966
        self.mat_clad.comp['260580'] = 1.28788E-05
        self.mat_clad = self.mat_clad*(clad_density/self.mat_clad.mass())

        self.enrichment = self.mat_absorber.comp[
            '50100']*10./(self.mat_absorber.comp['50100']*10. + self.mat_absorber.comp['50110']*11.)

        # define geometry
        self.absorber_dr = 0.3161
        self.clad_dr = 0.0559
        self.matrix_dr = 0

        self.setPF(0.502001622)

    def volume(self):
        import math

        self.absorber_r = self.absorber_dr
        self.clad_r = self.absorber_r + self.clad_dr

        self.A = math.pi*math.pow(self.clad_r, 2.)/self.PF
        self.pitch = math.pow(self.A/(math.tan(math.pi/6.)*6.), (1./2.))*2.

        self.matrix_r = math.pow(self.A/math.pi, .5)
        self.matrix_dr = self.matrix_r - self.clad_r

        self.absorber_da = math.pi*math.pow(self.absorber_r, 2.)
        self.clad_da = math.pi*math.pow(self.clad_r, 2.) - self.absorber_da
        self.matrix_da = math.pi * \
            math.pow(self.matrix_r, 2.) - self.clad_da - self.absorber_da

        self.absorber_vf = self.absorber_da/self.A
        self.clad_vf = self.clad_da/self.A
        self.matmrix_vf = self.matrix_da/self.A

    def setPF(self, PF):

        import math

        self.PF = PF
        self.volume()

    def setEnrichment(self, enrichment):

        import mocup

        mat = mocup.material()

        if enrichment > 1. or enrichment < 0.:
            poop

        AW = 1./(enrichment/10. + (1.-enrichment)/11.)
        b10 = AW/10.*enrichment
        b11 = AW/11.*(1-enrichment)

        mat.comp['50100'] = b10*4.
        mat.comp['50110'] = b11*4.
        mat.comp['60120'] = 1

        self.mat_absorber = mat*(self.absorber_density/mat.mass())

    def setPitch(self, pitch):
        import math

        self.volume()
        self.A = 12./2.*(pitch/2.)**2.*math.tan(math.pi/6.)
        self.PF = (self.absorber_da + self.clad_da)/self.A
        self.pitch = pitch


