#!/usr/bin/env python


class TRISO:
	
	def __init__(self):
		
		import mocup
		import math

		self.dr_kernel = .0175
		self.dr_buffer = .01
		self.dr_iPyC   = .0035
		self.dr_SiC    = .0035
		self.dr_oPyC   = .0035
		self.dr_matrix = 0.013573935

		self.mat_kernel = mocup.material()
		self.mat_kernel.comp = {'60120':(1.182344E-02/.6022),'80160':(3.547031E-02/.6022),'922350':(4.753724E-03/.6022),'922380':(1.889315E-02/.6022),'541350':1e-30}
		self.mat_buffer = mocup.material()
		self.mat_buffer.comp = {'60120':(5.01376E-02/.6022)}
		self.mat_iPyC = mocup.material()
		self.mat_iPyC.comp = {'60120':(9.37573E-02/.6022)}
		self.mat_SiC = mocup.material()
		self.mat_SiC.comp = {'60120':(4.81909E-02/.6022),'140280':(4.81909E-02/.6022)}
		self.mat_oPyC = mocup.material()
		self.mat_oPyC.comp = {'60120':(9.37573E-02/.6022)}
		self.mat_matrix = mocup.material()
		self.mat_matrix.comp = {'60120':(8.02201E-02/.6022)}

		self.max_k_err = .0001

	def volume(self):
		#This method finds the volumes of all the constituents of the TRISO particle
		import math
		
		self.r_kernel =                self.dr_kernel
		self.r_buffer = self.r_kernel + self.dr_buffer
		self.r_iPyC   = self.r_buffer + self.dr_iPyC
		self.r_SiC    = self.r_iPyC   + self.dr_SiC
		self.r_oPyC   = self.r_SiC    + self.dr_oPyC
		self.r_matrix = self.r_oPyC   + self.dr_matrix

		self.dv_kernel = 4./3.*math.pi*math.pow(self.r_kernel,3.)
		self.dv_buffer = 4./3.*math.pi*(math.pow(self.r_buffer,3.) - math.pow(self.r_kernel,3.))
		self.dv_iPyC   = 4./3.*math.pi*(math.pow(self.r_iPyC,3.)   - math.pow(self.r_buffer,3.))
		self.dv_SiC    = 4./3.*math.pi*(math.pow(self.r_SiC,3.)    - math.pow(self.r_iPyC,3.))
		self.dv_oPyC   = 4./3.*math.pi*(math.pow(self.r_oPyC,3.)   - math.pow(self.r_SiC,3.))
		self.dv_matrix = 4./3.*math.pi*(math.pow(self.r_matrix,3.) - math.pow(self.r_oPyC,3.))
		self.V         = self.dv_kernel + self.dv_buffer + self.dv_iPyC + self.dv_SiC + self.dv_oPyC + self.dv_matrix
		self.hpitch    = math.pow(self.V,(1./3.))/2.

	def PackingFraction(self):
		#This method updates TRISO fuel geometry

		import math

		self.volume()
		self.PF = math.pow((self.r_oPyC/self.r_matrix),3.)		

	def setPF(self, PF):
		#This method sets the TRISO fuel design to reflect a user defined packing fraction

		import math

		self.volume()
		
		self.PF = PF
		self.V = 4./3.*math.pi*math.pow(self.r_oPyC,3.)/PF
		self.r_matrix = math.pow((self.V/(4./3.*math.pi)),(1./3.))
		self.dr_matrix = self.r_matrix - self.r_oPyC
		self.dv_matrix = self.V - 4./3.*math.pi*math.pow(self.r_oPyC,3.)

	def homogenize(self):

		self.PackingFraction()
		self.mat_homogenized = self.mat_kernel*self.dv_kernel + self.mat_buffer*self.dv_buffer + self.mat_iPyC*self.dv_iPyC + self.mat_SiC*self.dv_SiC + self.mat_oPyC*self.dv_oPyC + self.mat_matrix*self.dv_matrix
		self.mat_homogenized = self.mat_homogenized*(1/self.V)

		return self.mat_homogenized

	def CHM(self):
		self.homogenize()
		return (self.mat_homogenized.comp['60120']/self.mat_homogenized.heavy_moles())

	def conductivities(self):
		
		import math

		self.homogenize()

		# conductivities W/mK
		self.k_kernel = 3.7
		self.k_buffer = 0.5
		self.k_iPyC   = 4.0
		self.k_SiC    = 16.0
		self.k_oPyC   = 4.0
		self.k_matrix = 15
		
		k = (self.k_kernel*self.dv_kernel + self.k_buffer*self.dv_buffer + self.k_iPyC*self.dv_iPyC + self.k_SiC*self.dv_SiC + self.k_oPyC*self.dv_oPyC + self.k_matrix*self.dv_matrix)/self.V

		ks = {}
		ks['min'] = 0.8*k
		ks['mid'] = 1.0*k
		ks['max'] = 1.2*k

		k_errs = {}
		k_errs['min'] = self.k_err(ks['min'])
		k_errs['max'] = self.k_err(ks['max'])
		k_errs['mid'] = 1.

		while math.fabs(k_errs['mid']) > self.max_k_err:
			ks['mid'] = (ks['min']-ks['max'])/(k_errs['min']-k_errs['max'])*(-k_errs['max'])+ks['max']
			k_errs['mid'] = self.k_err(ks['mid'])

			if ks['mid'] > ks['max']:
				ks['min'] = ks['max']
				k_errs['min'] = k_errs['max']
				ks['max'] = ks['mid']
				k_errs['max'] = k_errs['mid']
			elif ks['mid'] < ks['min']:
				ks['max'] = ks['min']
				k_errs['max'] = k_errs['min']
				ks['min'] = ks['mid']
				k_errs['min'] = k_errs['mid']
			elif k_errs['mid'] < 0.0:
				ks['max'] = ks['mid']
				k_errs['max'] = k_errs['mid']
			else:
				ks['min'] = ks['mid']
				k_errs['min'] = k_errs['mid']

		self.k_TRISO = ks['mid']

		return  ks['mid']
				
	def k_err(self,k):

		v = self.dv_kernel/self.V
		ki = self.k_kernel
		sum = v*(ki - k)/(ki + 2*k)
		v = self.dv_buffer/self.V
		ki = self.k_buffer
		sum += v*(ki - k)/(ki + 2*k)
		v = self.dv_iPyC/self.V	
		ki = self.k_iPyC
		sum += v*(ki - k)/(ki + 2*k)
		v = self.dv_SiC/self.V	
		ki = self.k_SiC
		sum += v*(ki - k)/(ki + 2*k)
		v = self.dv_oPyC/self.V
		ki = self.k_oPyC
		sum += v*(ki - k)/(ki + 2*k)
		v = self.dv_matrix/self.V
		ki = self.k_matrix
		sum += v*(ki - k)/(ki + 2*k)

		return sum

	def HT(self, Q):
		
		#This method returns the average temperature perturbation in the fuel kernel, ave_T and the maximum temperature perturbation in the fuel kernel, max_T

		import math

		# Q power per kernel (W)
		self.conductivities()
		
		# assume constant power density (W/m3) in kernel
		q = Q/(self.dv_kernel*1e-6)
		
		dt_kernel = q/(6*self.k_kernel)*math.pow((self.r_kernel*.01),2.)
		at_kernel = q/(10*self.k_kernel)*math.pow((self.r_kernel*.01),2.)

		r1 = 0.01*self.r_kernel
		r2 = 0.01*self.r_buffer
		k  = self.k_buffer

		dt_buffer = Q*(1./r1 - 1./r2)/(4*math.pi*k)
		at_buffer = Q*(3/(8*math.pi*k)*(math.pow(r2,2.) - math.pow(r1,2.))/(math.pow(r2,3.) - math.pow(r1,3.))-1/(4.*math.pi*k*r2))

		r1 = 0.01*self.r_buffer
		r2 = 0.01*self.r_iPyC
		k = self.k_iPyC
		
		dt_iPyC   = Q*(1./r1 - 1./r2)/(4*math.pi*k)
		at_iPyC   = Q*(3/(8*math.pi*k)*(math.pow(r2,2.) - math.pow(r1,2.))/(math.pow(r2,3.) - math.pow(r1,3.))-1/(4.*math.pi*k*r2))

		r1 = 0.01*self.r_iPyC
		r2 = 0.01*self.r_SiC
		k = self.k_SiC

		dt_SiC    = Q*(1./r1 - 1./r2)/(4*math.pi*k)
		at_SiC    = Q*(3/(8*math.pi*k)*(math.pow(r2,2.) - math.pow(r1,2.))/(math.pow(r2,3.) - math.pow(r1,3.))-1/(4.*math.pi*k*r2))

		r1 = 0.01*self.r_SiC
		r2 = 0.01*self.r_oPyC
		k = self.k_oPyC
		
		dt_oPyC   = Q*(1./r1 - 1./r2)/(4*math.pi*k)
		at_oPyC   = Q*(3/(8*math.pi*k)*(math.pow(r2,2.) - math.pow(r1,2.))/(math.pow(r2,3.) - math.pow(r1,3.))-1/(4.*math.pi*k*r2))

		r1 = 0.01*self.r_oPyC
		r2 = 0.01*self.r_matrix
		k = self.k_matrix
		
		dt_matrix = Q*(1./r1 - 1./r2)/(4*math.pi*k)
		at_matrix = Q*(3/(8*math.pi*k)*(math.pow(r2,2.) - math.pow(r1,2.))/(math.pow(r2,3.) - math.pow(r1,3.))-1/(4.*math.pi*k*r2))

		T_matrix = at_matrix
		T_oPyC   = dt_matrix + at_oPyC
		T_SiC    = dt_matrix + dt_oPyC + at_SiC
		T_iPyC   = dt_matrix + dt_oPyC + dt_SiC + at_iPyC
		T_buffer = dt_matrix + dt_oPyC + dt_SiC + dt_iPyC + at_buffer
		T_kernel = dt_matrix + dt_oPyC + dt_SiC + dt_iPyC + dt_buffer + at_kernel 

		offset = -(T_matrix*self.dv_matrix + T_oPyC*self.dv_oPyC + T_SiC*self.dv_SiC + T_iPyC*self.dv_iPyC + T_buffer*self.dv_buffer + T_kernel*self.dv_kernel)/self.V

		self.ave_T = T_kernel + offset
		self.max_T = T_kernel - at_kernel + dt_kernel + offset
		
		return self.ave_T, self.max_T

	def matrix(self):

		self.homogenize()
		self.mat_hmatrix = (self.mat_homogenized*self.V + self.mat_kernel*-self.dv_kernel)*(1./(self.V - self.dv_kernel))
		for nucl in list(self.mat_hmatrix.comp.keys()):
			if nucl in ['60120','60130','140280']:
				pass
			else:
				del self.mat_hmatrix.comp[nucl]

		return self.mat_hmatrix
