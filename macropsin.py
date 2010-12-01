import scipy.integrate as si


def smooth(t, start, stop):
    """
    Smoothly join (start, 0) to (stop, 1).
    """
    if t < start:
        t = 0
    else:
        if t > stop:
            t = 1
        else:
            t = (t - start)/(stop - start)
            
    return sin(t * pi/2)**2


    
class Constants:

	def __init__(self):
		self.MAGNETIZATION = 1.816 # T (Co @ 0 K)
		self.GAMMA_E = 176.0859770 # ns-1 T-1
		self.KB = 1.38e-23 # J K-1
		self.VOLUME = 1.13e-25 # m3
		self.MU0 = pi * 4e-7 # T m A-1
		self.NANOSECONDS = self.MAGNETIZATION * self.GAMMA_E # time unit  = 0.00312 ns
		self.GHZ = 1. / self.NANOSECONDS # freq. unit = 320 GHz
		self.TESLA = 1. / self.MAGNETIZATION # field unit = 1.816 T
		self.KELVIN = self.TESLA * self.KB
		self.H_K = 0.3 * self.TESLA
KONSTANTS = Constants()


class MacrospinParameters:

	def __init__(self, gamma=1, alpha=0.01, K=array([0., 0., 1.]), initial_magnetization=array([0., 0., 1.])):
		self.gamma = gamma
		self.alpha = alpha
		self.K = K * KONSTANTS.H_K/2.
		self.initial_magnetization = initial_magnetization * KONSTANTS.TESLA
		self.magnetization = copy(self.initial_magnetization)
		self.magnetization_trajectory = None
		self.energy = 0.

	def view(self):
		print("gamma: " + str(self.gamma))
		print("alpha: " + str(self.alpha))
		print("K: " + str(self.K))
		print("initial_magnetization: " + str(self.initial_magnetization))
		print("magnetization: " + str(self.magnetization))
		print("magnetization_trajectory: " + str(self.magnetization_trajectory))
		print("energy: " + str(self.energy))
		

	def normalize(self):
		norm = sqrt(sum(power(self.initial_magnetization, 2)))
		self.magnetization = divide(self.initial_magnetization, float(norm))


	def energy(self, h_dc, m):
		return -1. * dot(h_dc, m) - dot(self.K, power(m, 2)) #+ h_dc[1]



class MicrowaveParameters:

	def __init__(self, amplitude = 0.1, length = 1., rise = 0.1, phi = array([0., 3., 0.])):
		self.amplitude = amplitude * KONSTANTS.TESLA
		self.length = length * KONSTANTS.NANOSECONDS
		self.rise = rise * KONSTANTS.NANOSECONDS
		self.phi = phi
		self.env = 0.
		self.h_ac = array([0., 0., 0.])


	def view(self):
		print("amplitude: " + str(self.amplitude))
		print("length: " + str(self.length))
		print("rise: " + str(self.rise))
		print("phi: " + str(self.phi))
		print("phase: " + str(self.phase))
		print("env: " + str(self.env))
		print("h_ac: " + str(self.h_ac))
	

	def microwave_field(self, t):
		self.env = smooth(t, 0., self.rise) * (1 - smooth(t, self.length - self.rise, self.length))
		self.phase = dot(self.phi, power(t, range(len(self.phi))))
		self.h_ac[0] = self.amplitude * self.env * cos(self.phase*pi/180.)
		self.h_ac[1] = self.amplitude * self.env * sin(self.phase*pi/180.)
		self.h_ac[2] = 0.



class FieldParameters:

	def __init__(self, static_h=array([0., 0.1, +0.5])):
		self.static_h = static_h * KONSTANTS.H_K
		self.h_K =array([0., 0., 0.])
		self.h_dc = array([0., 0., 0.])
		self.h_st = array([0., 0., 0.])
		self.h_eff = array([0., 0., 0.])

	def view(self):
		print("static_h: " + str(self.static_h))
		print("h_K: " + str(self.h_K))
		print("h_dc: " + str(self.h_dc))
		print("h_st: " + str(self.h_st))
		print("h_eff: " + str(self.h_eff))

	def anisotropy_field(self, K, m):
		self.h_K = 2 * dot(K, m)


	def static_field(self, m, t):
		env = smooth(t, -10.* KONSTANTS.NANOSECONDS, -5.* KONSTANTS.NANOSECONDS)
		self.h_dc = env * self.static_h
		self.h_st = self.h_K + self.h_dc
		h_eff_m = sum(dot(self.h_st, m))
		self.h_st = self.h_st - h_eff_m * m


	def sumup_fields(self, h_ac):
		self.h_eff = self.h_st + h_ac




class MacrospinIntegration():

	def __init__(self, macrospin = MacrospinParameters(), field = FieldParameters(), rf = MicrowaveParameters(), length = 5., steptime = 0.05):

		self.macrospin = macrospin
		self.field = field
		self.rf = rf
		
		self.length = length * KONSTANTS.NANOSECONDS
		self.steptime = steptime * KONSTANTS.NANOSECONDS
		self.integration_length = self.length + 1. * KONSTANTS.NANOSECONDS
		self.time_sequence = arange(-10.* KONSTANTS.NANOSECONDS, self.integration_length, self.steptime)

		self.gamma = 0.
		self.alpha = 0.

		self.history = []
		
		self.t = []

		self.mx = []
		self.my = []
		self.mz = []

		self.hdcx = []
		self.hdcy = []
		self.hdcz = []

		self.hKx = []
		self.hKy = []
		self.hKz = []

		self.hacx = []
		self.hacy = []
		self.hacz = []

		self.heffx = []
		self.heffy = []
		self.heffz = []

		self.energy = []

	def view(self):
		print("macrospin: " + str(self.macrospin))
		print("field: " + str(self.field))
		print("rf: " + str(self.rf))
		print("length: " + str(self.length))
		print("steptime: " + str(self.steptime))
		print("integration_length: " + str(self.integration_length))
		print("time_sequence: " + str(self.time_sequence))
		print("gamma: " + str(self.gamma))
		print("alpha: " + str(self.alpha))

	def fast_relax(self, t):
		if t < 0:
			self.gamma = 1.e5
			self.alpha = 1.e4
		else:
			self.gamma = self.macrospin.gamma
			self.alpha = self.macrospin.alpha


	def cutoff_field(self, m, t):
		if m[2] < 0:
			self.integration_length = t

	
	def set_system_state(self, m, t):
		self.fast_relax(t)
		# self.cutoff_field(m, t)
		self.field.anisotropy_field(self.macrospin.K, m)
		self.field.static_field(m, t)
		self.rf.microwave_field(t)
		self.field.sumup_fields(self.rf.h_ac)
		self.t.append(t/KONSTANTS.NANOSECONDS)
		self.mx.append(m[0]/KONSTANTS.TESLA)
		self.my.append(m[1]/KONSTANTS.TESLA)
		self.mz.append(m[2]/KONSTANTS.TESLA)
		self.hdcx.append(self.field.h_dc[0]/KONSTANTS.TESLA)
		self.hdcy.append(self.field.h_dc[1]/KONSTANTS.TESLA)
		self.hdcz.append(self.field.h_dc[2]/KONSTANTS.TESLA)
		self.heffx.append(self.field.h_eff[0]/KONSTANTS.TESLA)
		self.heffy.append(self.field.h_eff[1]/KONSTANTS.TESLA)
		self.heffz.append(self.field.h_eff[2]/KONSTANTS.TESLA)
		self.hacx.append(self.rf.h_ac[0]/KONSTANTS.TESLA)
		self.hacy.append(self.rf.h_ac[1]/KONSTANTS.TESLA)
		self.hacz.append(self.rf.h_ac[2]/KONSTANTS.TESLA)
		# self.energy.append(self.macrospin.energy(self.field.h_dc, m))



#################################################################################

def LLG(m, t, integrate):
	integrate.set_system_state(m, t)
	return integrate.gamma * cross(m, integrate.field.h_eff - integrate.alpha * cross(m, cross(m, integrate.field.h_eff)))/(1 + integrate.alpha**2)
		

def do_integration(integrate):
	integrate.macrospin.normalize()
	integrate.macrospin.magnetization_trajectory = si.odeint(func=LLG, y0=integrate.macrospin.magnetization, t=integrate.time_sequence, args=(integrate, ))
