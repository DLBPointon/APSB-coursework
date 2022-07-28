import sbmltoodepy.modelclasses
from scipy.integrate import odeint
import numpy as np
import operator
import math

class LeloupModel(sbmltoodepy.modelclasses.Model):

	def __init__(self):

		self.p = {} #Dictionary of model parameters
		self.p['vsP'] = sbmltoodepy.modelclasses.Parameter(1.1, 'vsP', True, metadata = sbmltoodepy.modelclasses.SBMLMetadata("vsP"))
		self.p['vmP'] = sbmltoodepy.modelclasses.Parameter(1.0, 'vmP', True, metadata = sbmltoodepy.modelclasses.SBMLMetadata("vmP"))
		self.p['KmP'] = sbmltoodepy.modelclasses.Parameter(0.2, 'KmP', True, metadata = sbmltoodepy.modelclasses.SBMLMetadata("KmP"))
		self.p['KIP'] = sbmltoodepy.modelclasses.Parameter(1.0, 'KIP', True, metadata = sbmltoodepy.modelclasses.SBMLMetadata("KIP"))
		self.p['Pt'] = sbmltoodepy.modelclasses.Parameter(None, 'Pt', False, metadata = sbmltoodepy.modelclasses.SBMLMetadata("Pt"))
		self.p['ksP'] = sbmltoodepy.modelclasses.Parameter(0.9, 'ksP', True, metadata = sbmltoodepy.modelclasses.SBMLMetadata("ksP"))
		self.p['vdP'] = sbmltoodepy.modelclasses.Parameter(2.2, 'vdP', True, metadata = sbmltoodepy.modelclasses.SBMLMetadata("vdP"))
		self.p['KdP'] = sbmltoodepy.modelclasses.Parameter(0.2, 'KdP', True, metadata = sbmltoodepy.modelclasses.SBMLMetadata("KdP"))
		self.p['vsT'] = sbmltoodepy.modelclasses.Parameter(1.0, 'vsT', True, metadata = sbmltoodepy.modelclasses.SBMLMetadata("vsT"))
		self.p['vmT'] = sbmltoodepy.modelclasses.Parameter(0.7, 'vmT', True, metadata = sbmltoodepy.modelclasses.SBMLMetadata("vmT"))
		self.p['KmT'] = sbmltoodepy.modelclasses.Parameter(0.2, 'KmT', True, metadata = sbmltoodepy.modelclasses.SBMLMetadata("KmT"))
		self.p['KIT'] = sbmltoodepy.modelclasses.Parameter(1.0, 'KIT', True, metadata = sbmltoodepy.modelclasses.SBMLMetadata("KIT"))
		self.p['ksT'] = sbmltoodepy.modelclasses.Parameter(0.9, 'ksT', True, metadata = sbmltoodepy.modelclasses.SBMLMetadata("ksT"))
		self.p['vdT'] = sbmltoodepy.modelclasses.Parameter(3.0, 'vdT', True, metadata = sbmltoodepy.modelclasses.SBMLMetadata("vdT"))
		self.p['KdT'] = sbmltoodepy.modelclasses.Parameter(0.2, 'KdT', True, metadata = sbmltoodepy.modelclasses.SBMLMetadata("KdT"))
		self.p['kdC'] = sbmltoodepy.modelclasses.Parameter(0.01, 'kdC', True, metadata = sbmltoodepy.modelclasses.SBMLMetadata("kdC"))
		self.p['kdN'] = sbmltoodepy.modelclasses.Parameter(0.01, 'kdN', True, metadata = sbmltoodepy.modelclasses.SBMLMetadata("kdN"))
		self.p['k1'] = sbmltoodepy.modelclasses.Parameter(0.8, 'k1', True, metadata = sbmltoodepy.modelclasses.SBMLMetadata("k1"))
		self.p['k2'] = sbmltoodepy.modelclasses.Parameter(0.2, 'k2', True, metadata = sbmltoodepy.modelclasses.SBMLMetadata("k2"))
		self.p['k3'] = sbmltoodepy.modelclasses.Parameter(1.2, 'k3', True, metadata = sbmltoodepy.modelclasses.SBMLMetadata("k3"))
		self.p['k4'] = sbmltoodepy.modelclasses.Parameter(0.6, 'k4', True, metadata = sbmltoodepy.modelclasses.SBMLMetadata("k4"))
		self.p['kd'] = sbmltoodepy.modelclasses.Parameter(0.01, 'kd', True, metadata = sbmltoodepy.modelclasses.SBMLMetadata("kd"))
		self.p['V1P'] = sbmltoodepy.modelclasses.Parameter(8.0, 'V1P', True, metadata = sbmltoodepy.modelclasses.SBMLMetadata("V1P"))
		self.p['V1T'] = sbmltoodepy.modelclasses.Parameter(8.0, 'V1T', True, metadata = sbmltoodepy.modelclasses.SBMLMetadata("V1T"))
		self.p['V2P'] = sbmltoodepy.modelclasses.Parameter(1.0, 'V2P', True, metadata = sbmltoodepy.modelclasses.SBMLMetadata("V2P"))
		self.p['V2T'] = sbmltoodepy.modelclasses.Parameter(1.0, 'V2T', True, metadata = sbmltoodepy.modelclasses.SBMLMetadata("V2T"))
		self.p['V3P'] = sbmltoodepy.modelclasses.Parameter(8.0, 'V3P', True, metadata = sbmltoodepy.modelclasses.SBMLMetadata("V3P"))
		self.p['V3T'] = sbmltoodepy.modelclasses.Parameter(8.0, 'V3T', True, metadata = sbmltoodepy.modelclasses.SBMLMetadata("V3T"))
		self.p['V4P'] = sbmltoodepy.modelclasses.Parameter(1.0, 'V4P', True, metadata = sbmltoodepy.modelclasses.SBMLMetadata("V4P"))
		self.p['V4T'] = sbmltoodepy.modelclasses.Parameter(1.0, 'V4T', True, metadata = sbmltoodepy.modelclasses.SBMLMetadata("V4T"))
		self.p['K1P'] = sbmltoodepy.modelclasses.Parameter(2.0, 'K1P', True, metadata = sbmltoodepy.modelclasses.SBMLMetadata("K1P"))
		self.p['K1T'] = sbmltoodepy.modelclasses.Parameter(2.0, 'K1T', True, metadata = sbmltoodepy.modelclasses.SBMLMetadata("K1T"))
		self.p['K2P'] = sbmltoodepy.modelclasses.Parameter(2.0, 'K2P', True, metadata = sbmltoodepy.modelclasses.SBMLMetadata("K2P"))
		self.p['K2T'] = sbmltoodepy.modelclasses.Parameter(2.0, 'K2T', True, metadata = sbmltoodepy.modelclasses.SBMLMetadata("K2T"))
		self.p['K3P'] = sbmltoodepy.modelclasses.Parameter(2.0, 'K3P', True, metadata = sbmltoodepy.modelclasses.SBMLMetadata("K3P"))
		self.p['K3T'] = sbmltoodepy.modelclasses.Parameter(2.0, 'K3T', True, metadata = sbmltoodepy.modelclasses.SBMLMetadata("K3T"))
		self.p['K4P'] = sbmltoodepy.modelclasses.Parameter(2.0, 'K4P', True, metadata = sbmltoodepy.modelclasses.SBMLMetadata("K4P"))
		self.p['K4T'] = sbmltoodepy.modelclasses.Parameter(2.0, 'K4T', True, metadata = sbmltoodepy.modelclasses.SBMLMetadata("K4T"))
		self.p['n'] = sbmltoodepy.modelclasses.Parameter(4.0, 'n', True, metadata = sbmltoodepy.modelclasses.SBMLMetadata("n"))

		self.c = {} #Dictionary of compartments
		self.c['Compartment'] = sbmltoodepy.modelclasses.Compartment(1.0, 3, True, metadata = sbmltoodepy.modelclasses.SBMLMetadata(""))

		self.s = {} #Dictionary of chemical species
		self.s['MP'] = sbmltoodepy.modelclasses.Species(0.0614368, 'Concentration', self.c['Compartment'], False, constant = False, metadata = sbmltoodepy.modelclasses.SBMLMetadata("MP"))
		self.s['MP']._modifiedBy = 1
		self.s['CN'] = sbmltoodepy.modelclasses.Species(1.34728, 'Concentration', self.c['Compartment'], False, constant = False, metadata = sbmltoodepy.modelclasses.SBMLMetadata("CN"))
		self.s['CN']._modifiedBy = 10
		self.s['C'] = sbmltoodepy.modelclasses.Species(0.207614, 'Concentration', self.c['Compartment'], False, constant = False, metadata = sbmltoodepy.modelclasses.SBMLMetadata("C"))
		self.s['C']._modifiedBy = 9
		self.s['T2'] = sbmltoodepy.modelclasses.Species(0.0145428, 'Concentration', self.c['Compartment'], False, constant = False, metadata = sbmltoodepy.modelclasses.SBMLMetadata("T2"))
		self.s['T2']._modifiedBy = 8
		self.s['T1'] = sbmltoodepy.modelclasses.Species(0.0213384, 'Concentration', self.c['Compartment'], False, constant = False, metadata = sbmltoodepy.modelclasses.SBMLMetadata("T1"))
		self.s['T1']._modifiedBy = 7
		self.s['T0'] = sbmltoodepy.modelclasses.Species(0.0217261, 'Concentration', self.c['Compartment'], False, constant = False, metadata = sbmltoodepy.modelclasses.SBMLMetadata("T0"))
		self.s['T0']._modifiedBy = 6
		self.s['MT'] = sbmltoodepy.modelclasses.Species(0.0860342, 'Concentration', self.c['Compartment'], False, constant = False, metadata = sbmltoodepy.modelclasses.SBMLMetadata("MT"))
		self.s['MT']._modifiedBy = 5
		self.s['P0'] = sbmltoodepy.modelclasses.Species(0.0169928, 'Concentration', self.c['Compartment'], False, constant = False, metadata = sbmltoodepy.modelclasses.SBMLMetadata("P0"))
		self.s['P0']._modifiedBy = 2
		self.s['P1'] = sbmltoodepy.modelclasses.Species(0.0141356, 'Concentration', self.c['Compartment'], False, constant = False, metadata = sbmltoodepy.modelclasses.SBMLMetadata("P1"))
		self.s['P1']._modifiedBy = 3
		self.s['P2'] = sbmltoodepy.modelclasses.Species(0.0614368, 'Concentration', self.c['Compartment'], False, constant = False, metadata = sbmltoodepy.modelclasses.SBMLMetadata("P2"))
		self.s['P2']._modifiedBy = 4

		self.r = {} #Dictionary of reactions

		self.f = {} #Dictionary of function definitions
		self.time = 0

		self.AssignmentRules()



	def AssignmentRules(self):

		self.p['Pt'].value = self.s['P0'].concentration + self.s['P1'].concentration + self.s['P2'].concentration + self.s['C'].concentration + self.s['CN'].concentration

		return

	def RateMP(self):

		return self.p['vsP'].value * (self.p['KIP'].value**self.p['n'].value / (self.p['KIP'].value**self.p['n'].value + self.s['CN'].concentration**self.p['n'].value)) - (self.p['vmP'].value * (self.s['MP'].concentration / (self.p['KmP'].value + self.s['MP'].concentration)) + self.p['kd'].value * self.s['MP'].concentration)

	def RateP0(self):

		return self.p['ksP'].value * self.s['MP'].concentration + self.p['V2P'].value * (self.s['P1'].concentration / (self.p['K2P'].value + self.s['P1'].concentration)) - (self.p['V1P'].value * (self.s['P0'].concentration / (self.p['K1P'].value + self.s['P0'].concentration)) + self.p['kd'].value * self.s['P0'].concentration)

	def RateP1(self):

		return self.p['V1P'].value * (self.s['P0'].concentration / (self.p['K1P'].value + self.s['P0'].concentration)) + self.p['V4P'].value * (self.s['P2'].concentration / (self.p['K4P'].value + self.s['P2'].concentration)) - (self.p['V2P'].value * (self.s['P1'].concentration / (self.p['K2P'].value + self.s['P1'].concentration)) + self.p['V3P'].value * (self.s['P1'].concentration / (self.p['K3P'].value + self.s['P1'].concentration)) + self.p['kd'].value * self.s['P1'].concentration)

	def RateP2(self):

		return self.p['V3P'].value * (self.s['P1'].concentration / (self.p['K3P'].value + self.s['P1'].concentration)) + self.p['k4'].value * self.s['C'].concentration - (self.p['V4P'].value * (self.s['P2'].concentration / (self.p['K4P'].value + self.s['P2'].concentration)) + self.p['k3'].value * self.s['P2'].concentration * self.s['T2'].concentration + self.p['vdP'].value * (self.s['P2'].concentration / (self.p['KdP'].value + self.s['P2'].concentration)) + self.p['kd'].value * self.s['P2'].concentration)

	def RateMT(self):

		return self.p['vsT'].value * (self.p['KIT'].value**self.p['n'].value / (self.p['KIT'].value**self.p['n'].value + self.s['CN'].concentration**self.p['n'].value)) - (self.p['vmT'].value * (self.s['MT'].concentration / (self.p['KmT'].value + self.s['MT'].concentration)) + self.p['kd'].value * self.s['MT'].concentration)

	def RateT0(self):

		return self.p['ksT'].value * self.s['MT'].concentration + self.p['V2T'].value * (self.s['T1'].concentration / (self.p['K2T'].value + self.s['T1'].concentration)) - (self.p['V1T'].value * (self.s['T0'].concentration / (self.p['K1T'].value + self.s['T0'].concentration)) + self.p['kd'].value * self.s['T0'].concentration)

	def RateT1(self):

		return self.p['V1T'].value * (self.s['T0'].concentration / (self.p['K1T'].value + self.s['T0'].concentration)) + self.p['V4T'].value * (self.s['T2'].concentration / (self.p['K4T'].value + self.s['T2'].concentration)) - (self.p['V2T'].value * (self.s['T1'].concentration / (self.p['K2T'].value + self.s['T1'].concentration)) + self.p['V3T'].value * (self.s['T1'].concentration / (self.p['K3T'].value + self.s['T1'].concentration)) + self.p['kd'].value * self.s['T1'].concentration)

	def RateT2(self):

		return self.p['V3T'].value * (self.s['T1'].concentration / (self.p['K3T'].value + self.s['T1'].concentration)) + self.p['k4'].value * self.s['C'].concentration - (self.p['V4T'].value * (self.s['T2'].concentration / (self.p['K4T'].value + self.s['T2'].concentration)) + self.p['k3'].value * self.s['P2'].concentration * self.s['T2'].concentration + self.p['vdT'].value * (self.s['T2'].concentration / (self.p['KdT'].value + self.s['T2'].concentration)) + self.p['kd'].value * self.s['T2'].concentration)

	def RateC(self):

		return self.p['k3'].value * self.s['P2'].concentration * self.s['T2'].concentration + self.p['k2'].value * self.s['CN'].concentration - (self.p['k4'].value * self.s['C'].concentration + self.p['k1'].value * self.s['C'].concentration + self.p['kdC'].value * self.s['C'].concentration)

	def RateCN(self):

		return self.p['k1'].value * self.s['C'].concentration - (self.p['k2'].value * self.s['CN'].concentration + self.p['kdN'].value * self.s['CN'].concentration)

	def _SolveReactions(self, y, t):

		self.time = t
		self.s['MP'].amount, self.s['CN'].amount, self.s['C'].amount, self.s['T2'].amount, self.s['T1'].amount, self.s['T0'].amount, self.s['MT'].amount, self.s['P0'].amount, self.s['P1'].amount, self.s['P2'].amount = y
		self.AssignmentRules()

		rateRuleVector = np.array([ self.RateMP(), self.RateCN(), self.RateC(), self.RateT2(), self.RateT1(), self.RateT0(), self.RateMT(), self.RateP0(), self.RateP1(), self.RateP2()], dtype = np.float64)

		stoichiometricMatrix = np.array([[0.],[0.],[0.],[0.],[0.],[0.],[0.],[0.],[0.],[0.]], dtype = np.float64)

		reactionVelocities = np.array([0], dtype = np.float64)

		rateOfSpeciesChange = stoichiometricMatrix @ reactionVelocities + rateRuleVector

		return rateOfSpeciesChange

	def RunSimulation(self, deltaT, absoluteTolerance = 1e-12, relativeTolerance = 1e-6):

		finalTime = self.time + deltaT
		y0 = np.array([self.s['MP'].amount, self.s['CN'].amount, self.s['C'].amount, self.s['T2'].amount, self.s['T1'].amount, self.s['T0'].amount, self.s['MT'].amount, self.s['P0'].amount, self.s['P1'].amount, self.s['P2'].amount], dtype = np.float64)
		self.s['MP'].amount, self.s['CN'].amount, self.s['C'].amount, self.s['T2'].amount, self.s['T1'].amount, self.s['T0'].amount, self.s['MT'].amount, self.s['P0'].amount, self.s['P1'].amount, self.s['P2'].amount = odeint(self._SolveReactions, y0, [self.time, finalTime], atol = absoluteTolerance, rtol = relativeTolerance, mxstep=5000000)[-1]
		self.time = finalTime
		self.AssignmentRules()

