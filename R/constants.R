# Constants used in functions
# These are not exported
# Modified 5 Nov 2015 SDH

# Standard chemical formulas for mcomp
std.forms <- c(# Standard macromolecules
	       vfa = 'C2H4O2', protein = 'C5H7O2N', carbohydrate = 'C6H10O5', lipid = 'C57H104O6', 
	       lignin = 'C10H13O3', 
	       # Cell biomass 
	       biomass = 'C5H7O2N', 
	       # Organic acids
	       acetic = 'CH3COOH', lactic = 'C3H6O3', 
	       # Alcohols
	       ethanol = 'CH3CH2OH', 
	       cellulose = 'C6H10O5', 
	       # Carbohydrates
	       glucose = 'C6H12O6')

# Standard atomic weights
# These are from CIAAW (http://www.ciaaw.org/atomic-weights.htm), rounded to four digits if range was large 
atom.weights <- c(C = 12.01, H = 1.008, O = 16.00, N = 14.007, S = 32.06, P = 30.974, Na = 22.990, K = 39.098, Cl = 35.45)

# Molar gas volumes (mL/mol, at 273.15 K and 101.325 kPa (1 atm))
# From NIST (search webbook, then select fluid properties, then isobaric, etc.)
vol.mol <- c(CH4 = 22360.588, CO2 = 22263.009, N2 = 22403.863, H2 = 22427.978)
