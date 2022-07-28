import sbmltoodepy.modelclasses
from scipy.integrate import odeint
import numpy as np
import operator
import math

class BertozziModel(sbmltoodepy.modelclasses.Model):

	def __init__(self):

		self.p = {} #Dictionary of model parameters
		self.p['Pop_CA'] = sbmltoodepy.modelclasses.Parameter(39560000.0, 'Pop_CA', True, metadata = sbmltoodepy.modelclasses.SBMLMetadata("Pop_CA"))
		self.p['Ro_CA'] = sbmltoodepy.modelclasses.Parameter(2.7, 'Ro_CA', False, metadata = sbmltoodepy.modelclasses.SBMLMetadata("Ro_CA"))
		self.p['gamma_CA'] = sbmltoodepy.modelclasses.Parameter(0.14, 'gamma_CA', True, metadata = sbmltoodepy.modelclasses.SBMLMetadata("gamma_CA"))
		self.p['Io_CA'] = sbmltoodepy.modelclasses.Parameter(5.0, 'Io_CA', True, metadata = sbmltoodepy.modelclasses.SBMLMetadata("Io_CA"))
		self.p['Time'] = sbmltoodepy.modelclasses.Parameter(0.0, 'Time', False, metadata = sbmltoodepy.modelclasses.SBMLMetadata("Time"))
		self.p['Trigger_CA'] = sbmltoodepy.modelclasses.Parameter(1.0, 'Trigger_CA', True, metadata = sbmltoodepy.modelclasses.SBMLMetadata("Trigger_CA"))
		self.p['Trigger_NY'] = sbmltoodepy.modelclasses.Parameter(0.0, 'Trigger_NY', True, metadata = sbmltoodepy.modelclasses.SBMLMetadata("Trigger_NY"))
		self.p['Pop_NY'] = sbmltoodepy.modelclasses.Parameter(19540000.0, 'Pop_NY', True, metadata = sbmltoodepy.modelclasses.SBMLMetadata("Pop_NY"))
		self.p['Ro_NY'] = sbmltoodepy.modelclasses.Parameter(4.1, 'Ro_NY', False, metadata = sbmltoodepy.modelclasses.SBMLMetadata("Ro_NY"))
		self.p['gamma_NY'] = sbmltoodepy.modelclasses.Parameter(0.1, 'gamma_NY', True, metadata = sbmltoodepy.modelclasses.SBMLMetadata("gamma_NY"))
		self.p['Io_NY'] = sbmltoodepy.modelclasses.Parameter(29.0, 'Io_NY', True, metadata = sbmltoodepy.modelclasses.SBMLMetadata("Io_NY"))
		self.p['Pop'] = sbmltoodepy.modelclasses.Parameter(39560000.0, 'Pop', False, metadata = sbmltoodepy.modelclasses.SBMLMetadata("Pop"))
		self.p['Ro'] = sbmltoodepy.modelclasses.Parameter(2.7, 'Ro', False, metadata = sbmltoodepy.modelclasses.SBMLMetadata("Ro"))
		self.p['gamma'] = sbmltoodepy.modelclasses.Parameter(0.14, 'gamma', False, metadata = sbmltoodepy.modelclasses.SBMLMetadata("gamma"))
		self.p['Trigger_Lockdown'] = sbmltoodepy.modelclasses.Parameter(0.0, 'Trigger_Lockdown', True, metadata = sbmltoodepy.modelclasses.SBMLMetadata("Trigger_Lockdown"))
		self.p['Lockdown_CA_start'] = sbmltoodepy.modelclasses.Parameter(27.0, 'Lockdown_CA_start', True, metadata = sbmltoodepy.modelclasses.SBMLMetadata("Lockdown_CA_start"))
		self.p['Lockdown_CA_end'] = sbmltoodepy.modelclasses.Parameter(66.0, 'Lockdown_CA_end', True, metadata = sbmltoodepy.modelclasses.SBMLMetadata("Lockdown_CA_end"))
		self.p['Lockdown_NY_start'] = sbmltoodepy.modelclasses.Parameter(30.0, 'Lockdown_NY_start', True, metadata = sbmltoodepy.modelclasses.SBMLMetadata("Lockdown_NY_start"))
		self.p['Lockdown_NY_end'] = sbmltoodepy.modelclasses.Parameter(67.0, 'Lockdown_NY_end', True, metadata = sbmltoodepy.modelclasses.SBMLMetadata("Lockdown_NY_end"))
		self.p['Io'] = sbmltoodepy.modelclasses.Parameter(5.0, 'Io', False, metadata = sbmltoodepy.modelclasses.SBMLMetadata("Io"))
		self.p['Peak_Time'] = sbmltoodepy.modelclasses.Parameter(0.0, 'Peak_Time', False, metadata = sbmltoodepy.modelclasses.SBMLMetadata("Peak_Time"))
		self.p['ModelValue_3'] = sbmltoodepy.modelclasses.Parameter(5.0, 'ModelValue_3', True, metadata = sbmltoodepy.modelclasses.SBMLMetadata("Initial for Io_CA"))
		self.p['ModelValue_10'] = sbmltoodepy.modelclasses.Parameter(29.0, 'ModelValue_10', True, metadata = sbmltoodepy.modelclasses.SBMLMetadata("Initial for Io_NY"))
		self.p['ModelValue_16'] = sbmltoodepy.modelclasses.Parameter(66.0, 'ModelValue_16', True, metadata = sbmltoodepy.modelclasses.SBMLMetadata("Initial for Lockdown_CA_end"))
		self.p['ModelValue_15'] = sbmltoodepy.modelclasses.Parameter(27.0, 'ModelValue_15', True, metadata = sbmltoodepy.modelclasses.SBMLMetadata("Initial for Lockdown_CA_start"))
		self.p['ModelValue_18'] = sbmltoodepy.modelclasses.Parameter(67.0, 'ModelValue_18', True, metadata = sbmltoodepy.modelclasses.SBMLMetadata("Initial for Lockdown_NY_end"))
		self.p['ModelValue_17'] = sbmltoodepy.modelclasses.Parameter(30.0, 'ModelValue_17', True, metadata = sbmltoodepy.modelclasses.SBMLMetadata("Initial for Lockdown_NY_start"))
		self.p['ModelValue_0'] = sbmltoodepy.modelclasses.Parameter(39560000.0, 'ModelValue_0', True, metadata = sbmltoodepy.modelclasses.SBMLMetadata("Initial for Pop_CA"))
		self.p['ModelValue_7'] = sbmltoodepy.modelclasses.Parameter(19540000.0, 'ModelValue_7', True, metadata = sbmltoodepy.modelclasses.SBMLMetadata("Initial for Pop_NY"))
		self.p['ModelValue_1'] = sbmltoodepy.modelclasses.Parameter(2.7, 'ModelValue_1', True, metadata = sbmltoodepy.modelclasses.SBMLMetadata("Initial for Ro_CA"))
		self.p['ModelValue_8'] = sbmltoodepy.modelclasses.Parameter(4.1, 'ModelValue_8', True, metadata = sbmltoodepy.modelclasses.SBMLMetadata("Initial for Ro_NY"))
		self.p['ModelValue_5'] = sbmltoodepy.modelclasses.Parameter(1.0, 'ModelValue_5', True, metadata = sbmltoodepy.modelclasses.SBMLMetadata("Initial for Trigger_CA"))
		self.p['ModelValue_14'] = sbmltoodepy.modelclasses.Parameter(0.0, 'ModelValue_14', True, metadata = sbmltoodepy.modelclasses.SBMLMetadata("Initial for Trigger_Lockdown"))
		self.p['ModelValue_6'] = sbmltoodepy.modelclasses.Parameter(0.0, 'ModelValue_6', True, metadata = sbmltoodepy.modelclasses.SBMLMetadata("Initial for Trigger_NY"))
		self.p['ModelValue_2'] = sbmltoodepy.modelclasses.Parameter(0.14, 'ModelValue_2', True, metadata = sbmltoodepy.modelclasses.SBMLMetadata("Initial for gamma_CA"))
		self.p['ModelValue_9'] = sbmltoodepy.modelclasses.Parameter(0.1, 'ModelValue_9', True, metadata = sbmltoodepy.modelclasses.SBMLMetadata("Initial for gamma_NY"))

		self.c = {} #Dictionary of compartments
		self.c['USA___CA__NY'] = sbmltoodepy.modelclasses.Compartment(1.0, 3, True, metadata = sbmltoodepy.modelclasses.SBMLMetadata("USA - CA, NY"))

		self.s = {} #Dictionary of chemical species
		self.s['Infected'] = sbmltoodepy.modelclasses.Species(1.2639029322548e-07, 'Concentration', self.c['USA___CA__NY'], False, constant = False, metadata = sbmltoodepy.modelclasses.SBMLMetadata("Infected"))
		self.s['Recovered'] = sbmltoodepy.modelclasses.Species(0.0, 'Concentration', self.c['USA___CA__NY'], False, constant = False, metadata = sbmltoodepy.modelclasses.SBMLMetadata("Recovered"))
		self.s['Susceptible'] = sbmltoodepy.modelclasses.Species(0.999999999999997, 'Concentration', self.c['USA___CA__NY'], False, constant = False, metadata = sbmltoodepy.modelclasses.SBMLMetadata("Susceptible"))

		self.r = {} #Dictionary of reactions
		self.r['Susceptible_to_Infected'] = Susceptible_to_Infected(self)
		self.r['Infected_to_Recovered'] = Infected_to_Recovered(self)

		self.f = {} #Dictionary of function definitions
		self.f['rateOf'] = rateOf(self)
		self.f['Rate_Law_for_Susceptible_to_Infected'] = Rate_Law_for_Susceptible_to_Infected(self)
		self.time = 0

		self.AssignmentRules()



	def AssignmentRules(self):

		if self.time <= 0 :
			isConstantValue = self.p['ModelValue_3']._constant
			self.p['ModelValue_3']._constant = False
			self.p['ModelValue_3'].value = self.p['Io_CA'].value
			self.p['ModelValue_3']._constant = isConstantValue

		if self.time <= 0 :
			isConstantValue = self.p['ModelValue_10']._constant
			self.p['ModelValue_10']._constant = False
			self.p['ModelValue_10'].value = self.p['Io_NY'].value
			self.p['ModelValue_10']._constant = isConstantValue

		if self.time <= 0 :
			isConstantValue = self.p['ModelValue_16']._constant
			self.p['ModelValue_16']._constant = False
			self.p['ModelValue_16'].value = self.p['Lockdown_CA_end'].value
			self.p['ModelValue_16']._constant = isConstantValue

		if self.time <= 0 :
			isConstantValue = self.p['ModelValue_15']._constant
			self.p['ModelValue_15']._constant = False
			self.p['ModelValue_15'].value = self.p['Lockdown_CA_start'].value
			self.p['ModelValue_15']._constant = isConstantValue

		if self.time <= 0 :
			isConstantValue = self.p['ModelValue_18']._constant
			self.p['ModelValue_18']._constant = False
			self.p['ModelValue_18'].value = self.p['Lockdown_NY_end'].value
			self.p['ModelValue_18']._constant = isConstantValue

		if self.time <= 0 :
			isConstantValue = self.p['ModelValue_17']._constant
			self.p['ModelValue_17']._constant = False
			self.p['ModelValue_17'].value = self.p['Lockdown_NY_start'].value
			self.p['ModelValue_17']._constant = isConstantValue

		if self.time <= 0 :
			isConstantValue = self.p['ModelValue_0']._constant
			self.p['ModelValue_0']._constant = False
			self.p['ModelValue_0'].value = self.p['Pop_CA'].value
			self.p['ModelValue_0']._constant = isConstantValue

		if self.time <= 0 :
			isConstantValue = self.p['ModelValue_7']._constant
			self.p['ModelValue_7']._constant = False
			self.p['ModelValue_7'].value = self.p['Pop_NY'].value
			self.p['ModelValue_7']._constant = isConstantValue

		if self.time <= 0 :
			isConstantValue = self.p['ModelValue_1']._constant
			self.p['ModelValue_1']._constant = False
			self.p['ModelValue_1'].value = self.p['Ro_CA'].value
			self.p['ModelValue_1']._constant = isConstantValue

		if self.time <= 0 :
			isConstantValue = self.p['ModelValue_8']._constant
			self.p['ModelValue_8']._constant = False
			self.p['ModelValue_8'].value = self.p['Ro_NY'].value
			self.p['ModelValue_8']._constant = isConstantValue

		if self.time <= 0 :
			isConstantValue = self.p['ModelValue_5']._constant
			self.p['ModelValue_5']._constant = False
			self.p['ModelValue_5'].value = self.p['Trigger_CA'].value
			self.p['ModelValue_5']._constant = isConstantValue

		if self.time <= 0 :
			isConstantValue = self.p['ModelValue_14']._constant
			self.p['ModelValue_14']._constant = False
			self.p['ModelValue_14'].value = self.p['Trigger_Lockdown'].value
			self.p['ModelValue_14']._constant = isConstantValue

		if self.time <= 0 :
			isConstantValue = self.p['ModelValue_6']._constant
			self.p['ModelValue_6']._constant = False
			self.p['ModelValue_6'].value = self.p['Trigger_NY'].value
			self.p['ModelValue_6']._constant = isConstantValue

		self.p['Pop'].value = self.p['ModelValue_5'].value * self.p['ModelValue_0'].value + self.p['ModelValue_6'].value * self.p['ModelValue_7'].value

		self.p['Ro'].value = self.p['ModelValue_5'].value * self.p['Ro_CA'].value + self.p['ModelValue_6'].value * self.p['Ro_NY'].value

		self.p['Io'].value = self.p['ModelValue_5'].value * self.p['ModelValue_3'].value + self.p['ModelValue_6'].value * self.p['ModelValue_10'].value

		if self.time <= 0 :
			isConstantValue = self.s['Infected']._constant
			self.s['Infected']._constant = False
			self.s['Infected'].concentration = self.p['Io'].value / self.p['Pop'].value
			self.s['Infected']._constant = isConstantValue

		if self.time <= 0 :
			isConstantValue = self.s['Susceptible']._constant
			self.s['Susceptible']._constant = False
			self.s['Susceptible'].concentration = (self.p['Pop'].value - self.s['Infected'].concentration) / self.p['Pop'].value
			self.s['Susceptible']._constant = isConstantValue

		if self.time <= 0 :
			isConstantValue = self.p['ModelValue_2']._constant
			self.p['ModelValue_2']._constant = False
			self.p['ModelValue_2'].value = self.p['gamma_CA'].value
			self.p['ModelValue_2']._constant = isConstantValue

		self.p['Time'].value = self.p['ModelValue_2'].value * self.time

		if self.time <= 0 :
			isConstantValue = self.p['ModelValue_9']._constant
			self.p['ModelValue_9']._constant = False
			self.p['ModelValue_9'].value = self.p['gamma_NY'].value
			self.p['ModelValue_9']._constant = isConstantValue

		self.p['gamma'].value = self.p['ModelValue_5'].value * self.p['ModelValue_2'].value + self.p['ModelValue_6'].value * self.p['ModelValue_9'].value

		return

	def _SolveReactions(self, y, t):

		self.time = t
		self.s['Infected'].amount, self.s['Recovered'].amount, self.s['Susceptible'].amount = y
		self.AssignmentRules()

		rateRuleVector = np.array([ 0, 0, 0], dtype = np.float64)

		stoichiometricMatrix = np.array([[ 1,-1.],[ 0,1.],[-1,0.]], dtype = np.float64)

		reactionVelocities = np.array([self.r['Susceptible_to_Infected'](), self.r['Infected_to_Recovered']()], dtype = np.float64)

		rateOfSpeciesChange = stoichiometricMatrix @ reactionVelocities + rateRuleVector

		return rateOfSpeciesChange

	def RunSimulation(self, deltaT, absoluteTolerance = 1e-12, relativeTolerance = 1e-6):

		finalTime = self.time + deltaT
		y0 = np.array([self.s['Infected'].amount, self.s['Recovered'].amount, self.s['Susceptible'].amount], dtype = np.float64)
		self.s['Infected'].amount, self.s['Recovered'].amount, self.s['Susceptible'].amount = odeint(self._SolveReactions, y0, [self.time, finalTime], atol = absoluteTolerance, rtol = relativeTolerance, mxstep=5000000)[-1]
		self.time = finalTime
		self.AssignmentRules()

class Susceptible_to_Infected:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		if metadata:
			self.metadata = metadata
		else:
			self.metadata = sbmltoodepy.modelclasses.SBMLMetadata("Susceptible_to_Infected")

	def __call__(self):
		return self.parent.c['USA___CA__NY'].size * self.parent.f['Rate_Law_for_Susceptible_to_Infected'](self.parent.p['Ro'].value, self.parent.s['Infected'].concentration, self.parent.s['Susceptible'].concentration, self.parent.p['gamma'].value)

class Infected_to_Recovered:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		if metadata:
			self.metadata = metadata
		else:
			self.metadata = sbmltoodepy.modelclasses.SBMLMetadata("Infected_to_Recovered")

	def __call__(self):
		return self.parent.c['USA___CA__NY'].size * self.parent.p['gamma'].value * self.parent.s['Infected'].concentration

class rateOf:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		if metadata:
			self.metadata = metadata
		else:
			self.metadata = sbmltoodepy.modelclasses.SBMLMetadata("")
	def __call__(self, a):
		return NaN

class Rate_Law_for_Susceptible_to_Infected:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		if metadata:
			self.metadata = metadata
		else:
			self.metadata = sbmltoodepy.modelclasses.SBMLMetadata("Rate Law for Susceptible_to_Infected")
	def __call__(self, Ro, I, S, gamma):
		return gamma * Ro * I * S

