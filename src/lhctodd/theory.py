from .tools import __data_path__
from scipy import interpolate
import numpy as np

class width:
    # quark masses from PDG 2020
    _q_mass = [
        0.00216, # up
        0.00467, # down
        0.093,   # strange
        1.27,    # charm
        4.18,    # bottom
        172.76   # top
    ]
    # lepton masses: e mu tau from PDG 2020
    _l_mass = [0.000511, 0.105658, 1.77682]
    _vev = 246
    
    # fermion-Yukawa-Coupling function 
    _fyc = lambda mf: np.sqrt(2)*mf/246.0
    
    # fetch alpha_s values taken from NNPDF31
    _as_data = np.genfromtxt(
        str( __data_path__ / f"alpha_s.csv"),
        delimiter=","
    )
    _as = interpolate.interp1d(
        _as_data[:,0], # scale
        _as_data[:,1], # alpha_s value
        fill_value="extrapolate"
    )

    # Vector mediators
    @classmethod
    def vector_qq(cls, med_mass, g=1.0):
        """ Width of vector mediator decaying to quarks        
        """
        width = 0
        for mass in cls._q_mass:
            z = np.divide(mass, med_mass)**2
            z = z.astype(np.complex128)
            width += np.where(
                med_mass >= 2 * mass,
                g**2 * med_mass * np.sqrt(1 - 4*z) * (1 + 2*z) / (4*np.pi),
                0.0 
            )
        return np.abs(width)

    @classmethod
    def vector_ll(cls, med_mass, g=1.0):
        """Width of vector mediator decaying to leptons
        """
        width = 0
        for mass in cls._l_mass:
            z = np.divide(mass, med_mass)**2
            z = z.astype(np.complex128)
            width += np.where(
                med_mass >= 2 * mass,
                g**2 * med_mass * np.sqrt(1 - 4*z) * (1 + 2*z) / (12*np.pi),
                0.0 
            )
        return np.abs(width)

    @classmethod
    def vector_nn(cls, med_mass, g=1.0):
        """Width of vector mediator decaying to neutrinos
        """
        return g**2 * med_mass / (24*np.pi)

    @classmethod
    def vector_dm(cls, med_mass, chi_mass=1.0, g=1.0):
        """Width of vector mediator decaying to dark matter candidates 
        """
        z = np.divide(chi_mass, med_mass)**2
        z = z.astype(np.complex128)
        return np.abs(np.where(
            med_mass >= 2 * chi_mass, 
            g**2 * med_mass * np.sqrt(1 - 4*z) * (1 + 2*z) / (12*np.pi),
            0.0
        ))

    @classmethod
    def vector_total(cls, med_mass, chi_mass, g_q=0.25, g_chi=1.0, g_l=0.0):
        """Total width of the vector mediator
        """
        total  = cls.vector_qq(med_mass, g_q)
        total += cls.vector_ll(med_mass, g_l)
        total += cls.vector_nn(med_mass, g_l)
        total += cls.vector_dm(med_mass, chi_mass, g_chi)
        return total
    
    #Axial-Vector Mediators
    @classmethod
    def axial_qq(cls, med_mass, g = 1.0):
        width = 0
        for mass in cls._q_mass:
            z = np.divide(mass, med_mass)**2
            z = z.astype(np.complex128)
            width += np.where(
                med_mass >= 2 * mass,
                g**2 * med_mass * np.sqrt(1 - 4*z)**3 / (4*np.pi),
                0.0 
            )
        return np.abs(width)

    @classmethod
    def axial_ll(cls, med_mass, g = 1.0):
        width = 0
        for mass in cls._l_mass:
            z = np.divide(mass, med_mass)**2
            z = z.astype(np.complex128)
            width += np.where(
                med_mass >= 2 * mass,
                g**2 * med_mass * np.sqrt(1 - 4*z)**3 / (12*np.pi),
                0.0 
            )
        return np.abs(width)

    @classmethod
    def axial_nn(cls, med_mass, g = 1.0):
        return g**2 * med_mass / (24*np.pi)

    @classmethod
    def axial_dm(cls, med_mass, chi_mass=1.0, g = 1.0):
        z = np.divide(chi_mass, med_mass)**2
        z = z.astype(np.complex128)
        return np.abs(np.where(
            med_mass >= 2 * chi_mass, 
            g**2 * med_mass * np.sqrt(1 - 4*z)**3 / (12*np.pi),
            0.0
        ))

    @classmethod
    def axial_total(cls, med_mass, chi_mass, g_q=0.25, g_dm=1.0, g_l=0.0):
        total  = cls.axial_qq(med_mass, g_q)
        total += cls.axial_ll(med_mass, g_l)
        total += cls.axial_nn(med_mass, g_l)
        total += cls.axial_dm(med_mass, chi_mass, g_dm)
        return total

    # Scalar Mediators
    def form_factor_s(tau):
        tau = tau.astype(np.complex128)
        return tau * (
            1 + (1 - tau)* np.arctan(np.divide(1, np.sqrt(tau - 1)))**2
        )
    def from_factor_ps(tau):
        tau = tau.astype(np.complex128)
        return np.abs(tau * np.arctan(np.divide(1, np.sqrt(tau - 1)))**2)

    @classmethod
    def scalar_qq(cls, med_mass, g=1.0):
        width = 0
        for mass in cls._q_mass:
            z = np.divide(mass, med_mass)**2
            z = z.astype(np.complex128)
            width += np.where(
                med_mass >= 2 * mass,
                3 * g**2 * cls._fyc(mass)**2 * med_mass * np.sqrt(1 - 4*z)**3 / (16*np.pi),
                0.0
            )
        return np.abs(width)
    
    @classmethod
    def scalar_gg(cls, med_mass, g=1.0):
        z = np.divide(cls._q_mass[5], med_mass)**2
        z = z.astype(np.complex128)
        return np.abs(g**2 * med_mass**3 * cls._as(med_mass)**2 * np.where(
            med_mass >= 2 * cls._q_mass[5], 
            cls.form_factor_s(4*z)**2 / (32*np.pi**3 * cls._vev**2),
            0.0
        ))

    @classmethod
    def scalar_qq(cls, med_mass, g=1.0):
        width = 0
        for mass in cls._q_mass:
            z = np.divide(mass, med_mass)**2
            z = z.astype(np.complex128)
            width += np.where(
                med_mass >= 2 * mass, 
                3 * g**2 * cls._fyc(mass)**2 * med_mass * np.sqrt(1 - 4*z**2)**3 / (16*np.pi),
                0.0
            )
        return np.abs(width)

    @classmethod
    def scalar_dm(cls, med_mass, chi_mass=1.0, g=1.0):
        z = np.divide(chi_mass, med_mass)**2
        z = z.astype(np.complex128)
        return np.abs(
            g**2 * med_mass * np.where(
                med_mass >= 2 * chi_mass,
                np.sqrt(1 - 4 * z**2)**3 / (8 * np.pi),
                0.0
            )
        )
    
    @classmethod
    def scalar_total(cls, med_mass, chi_mass=1.0, g_q=0.25, g_dm=1.0):
        total = cls.scalar_dm(med_mass, chi_mass, g_dm)
        total += cls.scalar_gg(med_mass, g_q)
        total += cls.scalar_qq(med_mass, g_q)
        return total

    # pseudo-scalar mediators
    @classmethod
    def pseudo_scalar_gg(cls, med_mass, g=1.0):
        z = np.divide(cls._q_mass[5], med_mass)**2
        z = z.astype(np.complex128)
        return np.abs(g**2 * med_mass**3 * cls._as(med_mass)**2 * np.where(
            med_mass >= 2 * cls._q_mass[5], 
            cls.form_factor_ps(4*z)**2 / (32*np.pi**3 * cls._vev**2),
            0.0
        ))

    @classmethod
    def pseudo_scalar_qq(cls, med_mass, g=1.0):
        width = 0
        for mass in cls._q_mass:
            z = np.divide(mass, med_mass)**2
            z = z.astype(np.complex128)
            width += np.where(
                med_mass >= 2 * mass, 
                3 * g**2 * cls._fyc(mass)**2 * med_mass * np.sqrt(1 - 4*z**2)/ (16*np.pi),
                0.0
            )
        return np.abs(width)

    @classmethod
    def pseudo_scalar_dm(cls, med_mass, chi_mass=1.0, g=1.0):
        z = np.divide(chi_mass, med_mass)**2
        z = z.astype(np.complex128)
        return np.abs(
            g**2 * med_mass * np.where(
                med_mass >= 2 * chi_mass,
                np.sqrt(1 - 4 * z**2) / (8 * np.pi),
                0.0
            )
        )
    
    @classmethod
    def pseudo_scalar_total(cls, med_mass, chi_mass=1.0, g_q=0.25, g_g=0.0, g_dm=1.0):
        total = cls.pseudo_scalar_dm(med_mass, chi_mass, g_dm)
        total += cls.pseudo_scalar_gg(med_mass, g_q)
        total += cls.pseudo_scalar_qq(med_mass, g_q)
        return total

