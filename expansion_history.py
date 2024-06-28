from scipy.integrate import quad
import numpy as np
import argparse
from copy import deepcopy
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('Agg')


#######################################
# Create command line argument parser
#######################################

def create_parser():

        # Handle user input with argparse
    parser = argparse.ArgumentParser(
        description="Compute the expansion history.")

    parser.add_argument('--Omega_m',
        dest='Omega_m',
        default=0.3111,
        type=float,
        help='Matter overdensity Omega_m.')

    parser.add_argument('--Omega_r',
        dest='Omega_r',
#        default=4.203069657615112e-05,
        default=0,
        type=float,
        help='Radiation overdensity Omega_r h**2.')

    parser.add_argument('--Omega_k',
        dest='Omega_k',
        default=0,
        type=float,
        help='Curvature overdensity Omega_k.')

    parser.add_argument('--Omega_de',
        dest='Omega_de',
        default=None,
        type=float,
        help='Dark energy overdensity Omega_de.')

    parser.add_argument('--w0',
        dest='w0',
        default=-1,
        type=float,
        help='Dark energy equation of state at z=0.')

    parser.add_argument('--wa',
        dest='wa',
        default=0,
        type=float,
        help='Dark energy equation of state parameter defined by w(a) = w0 + wa(1-a).')

    parser.add_argument('--H0',
        dest='H0',
        default=67.66,
        type=float,
        help='Hubble parameter in km/s/Mpc.')

    parser.add_argument('--name',
        dest='name',
        default='User Specified',
        type=str,
        help='Cosmology name.')

    parser.add_argument('-v', '--verbose',
        dest='verbose',
        action='store_true',
        help='Print helpful information to the screen? (default: False)',
        default=False)

    parser.add_argument('--Concordance',
        dest='Concordance',
        action='store_true',
        help='Use Concordance Cosmology? (default: False)',
        default=False)

    parser.add_argument('--DESI+CMB+Pantheon',
        dest='DESI_CMB_Pantheon',
        action='store_true',
        help='Use DESI+CMB+Pantheon Cosmology? (default: False)',
        default=False)

    parser.add_argument('--DESI+CMB+Union3',
        dest='DESI_CMB_Union3',
        action='store_true',
        help='Use DESI+CMB+Union3 Cosmology? (default: False)',
        default=False)

    parser.add_argument('--DESI+CMB+DESY5',
        dest='DESI_CMB_DESY5',
        action='store_true',
        help='Use DESI+CMB+DESY5 Cosmology? (default: False)',
        default=False)

    return parser

class CosmologicalParameters:
	name = 'User Specified'
	H0 = 0
	Omega_m = 0
	Omega_r = 0
	Omega_de = 0
	Omega_k = 0
	w0 = 0
	wa = 0
	DESI_CMB_DESY5 = False
	DESI_CMB_Union3 = False
	DESI_CMB_Pantheon = False
	Concordance = False
	def __init__(self):
		self.name = 'User Specified'
		self.H0 = 0
		self.Omega_m = 0
		self.Omega_r = 0
		self.Omega_de = 0
		self.Omega_k = 0
		self.w0 = 0
		self.wa = 0
		self.DESI_CMB_DESY5 = False
		self.DESI_CMB_Union3 = False
		self.DESI_CMB_Pantheon = False
		self.Concordance = False
	def copy(self):
		return deepcopy(self)


def plot_scale_factor_vs_look_back_time(tlb,a,args=None):

	#get an LCDM baseline
	z = 1./a - 1.
	parser = create_parser()
	argslcdm = CosmologicalParameters()
	argsdcp  = argslcdm.copy()#parser.parse_args()
	argsdcu  = argslcdm.copy()
	argsdcd  = argslcdm.copy()


	#concordance
	argslcdm.Concordance = True
	argslcdm = special_cosmologies(argslcdm)
	argslcdm = adjust_omegas(argslcdm)
	tllcdm = t(z,args=argslcdm)
	H0_inv = 1./H0_in_inv_gyr(args=argslcdm)
	tllcdm/=H0_inv #H0 t


	#DESI+CMB+Pantheon
	argsdcp.DESI_CMB_Pantheon = True
	argsdcp = special_cosmologies(argsdcp)
	argsdcp = adjust_omegas(argsdcp)
	tldcp = t(z,args=argsdcp)
	H0_inv = 1./H0_in_inv_gyr(args=argsdcp)
	tldcp/=H0_inv #H0 t

	#DESI+CMB+Union3
	argsdcu.DESI_CMB_Union3 = True
	argsdcu = special_cosmologies(argsdcu)
	argsdcu = adjust_omegas(argsdcu)
	tldcu = t(z,args=argsdcu)
	H0_inv = 1./H0_in_inv_gyr(args=argsdcu)
	tldcu/=H0_inv #H0 t

	#DESI+CMB+DESY5
	argsdcd.DESI_CMB_DESY5 = True
	argsdcd = special_cosmologies(argsdcd)
	argsdcd = adjust_omegas(argsdcd)
	tldcd = t(z,args=argsdcd)
	H0_inv = 1./H0_in_inv_gyr(args=argsdcd)
	tldcd/=H0_inv #H0 t

	print_cosmology(argslcdm)
	print_cosmology(argsdcp)
	print_cosmology(argsdcu)
	print_cosmology(argsdcd)





	f,ax = plt.subplots(1,2,figsize=(14,7))

	plt.style.use('robertsons_rules')

	#expansion history
	ax[0].plot(tldcp,a,label=fr'{argsdcp.name}: $w_0$={argsdcp.w0:4.3f}, $w_a$={argsdcp.wa:4.3f}')
	ax[0].plot(tldcu,a,label=fr'{argsdcu.name}: $w_0$={argsdcu.w0:4.3f}, $w_a$={argsdcu.wa:4.3f}')
	ax[0].plot(tldcd,a,label=fr'{argsdcd.name}: $w_0$={argsdcd.w0:4.3f}, $w_a$={argsdcd.wa:4.3f}')
	ax[0].plot(tllcdm,a,label='LCDM',linestyle=':')
	ax[0].set_xlim([1.,0])
	ax[0].set_ylim([0,1])
	ax[0].set_xlabel(r'Look-back time $H_0 t$')
	ax[0].set_ylabel(r'Scale factor $a$')
	ax[0].legend(frameon=False)

	a_tldcp_lcdm = np.interp(tldcp,tllcdm,a)
	a_tldcu_lcdm = np.interp(tldcu,tllcdm,a)
	a_tldcd_lcdm = np.interp(tldcd,tllcdm,a)

	ax[1].plot(tldcp,(a/a_tldcp_lcdm),label=fr'{argsdcp.name}: $w_0$={argsdcp.w0:4.3f}, $w_a$={argsdcp.wa:4.3f}')
	ax[1].plot(tldcu,(a/a_tldcu_lcdm),label=fr'{argsdcu.name}: $w_0$={argsdcu.w0:4.3f}, $w_a$={argsdcu.wa:4.3f}')
	ax[1].plot(tldcd,(a/a_tldcd_lcdm),label=fr'{argsdcd.name}: $w_0$={argsdcd.w0:4.3f}, $w_a$={argsdcd.wa:4.3f}')
	ax[1].plot(tllcdm,(a/a),label='LCDM',linestyle=':')
	ax[1].set_xlim([0.8,0])
	ax[1].set_ylim([0.90,1.01])
	ax[1].set_xlabel(r'Look-back time $H_0 t$')
	ax[1].set_ylabel(r'Scale factor $a$ / Scale factor $a_\Lambda(H_0t)$')
	ax[1].legend(frameon=False)

	plt.savefig('expansion_history.png',bbox_inches='tight',dpi=300)

def adjust_omega_r(args):
	#convert between Omega_r h**2 and Omega_r
	#This is based on the Planck 2018 result
	#that z_eq = 3387 and Omega_m h^2 = 0.1420
	#so the default Omega_r = Omega_m/(1+z_eq).~4.203069657615112e-05
	args.Omega_r /= (args.H0/100)**2

	return args

def adjust_omega_de(args):
	if(args.Omega_de is None):
		#assume flat cosmology
		args.Omega_de = 1.0 - (args.Omega_m+args.Omega_k+args.Omega_r)
	return args

def adjust_omegas(argsin):
	args = adjust_omega_r(argsin)
	args = adjust_omega_de(args)
	return args


def print_cosmology(args):
 
	#print parameter values
	print(f'Cosmological Model: {args.name}')
	print(f'Hubble      parameter   H0       {args.H0}')
	print(f'Matter      overdensity Omega_m  {args.Omega_m}')
	print(f'Dark Energy overdensity Omega_de {args.Omega_de}')
	print(f'Radiation   overdensity Omega_r  {args.Omega_r}')
	print(f'Curvature   overdensity Omega_k  {args.Omega_k}')
	print(f'DE equation of state    w0       {args.w0}')
	print(f'DE equation of state    wa       {args.wa}')


def special_cosmologies(argsin):

	#see table 3 of https://arxiv.org/pdf/2404.03002

	if(argsin.Concordance):
		argsin.name = 'Concordance'
		argsin.Omega_m = 0.3
		argsin.H0 = 70.
		#argsin.Omega_r = 4.203069657615112e-05
		argsin.Omega_r = 0
		argsin.Omega_de = None #set later
		argsin.w0 = -1
		argsin.wa = 0
	if(argsin.DESI_CMB_Pantheon):
		argsin.name = 'DESI+CMB+Pantheon'
		argsin.Omega_m = 0.3085
		argsin.H0 = 68.03
		argsin.Omega_r = 0#4.203069657615112e-05
		argsin.Omega_de = None #set later
		argsin.w0 = -0.827
		argsin.wa = -0.75
	if(argsin.DESI_CMB_Union3):
		argsin.name = 'DESI+CMB+Union3'
		argsin.Omega_m = 0.3230
		argsin.H0 = 66.53
		argsin.Omega_r = 4.203069657615112e-05
		argsin.Omega_de = None #set later
		argsin.w0 = -0.65
		argsin.wa = -1.27
	if(argsin.DESI_CMB_DESY5):
		argsin.name = 'DESI+CMB+DESY5'
		argsin.Omega_m = 0.3230
		argsin.H0 = 67.24
		argsin.Omega_r = 4.203069657615112e-05
		argsin.Omega_de = None #set later
		argsin.w0 = -0.727
		argsin.wa = -1.05

	return argsin



def print_reference_histories():

	#see table 3 of https://arxiv.org/pdf/2404.03002

	w0       = -1.0
	wa       = 0.0
	print(f'LCDM has w0 = {w0}, wa = {wa}.\n')

	print(f'Here is the constrain table from DESI, see https://arxiv.org/pdf/2404.03002 :')
	print(f'Omega_m\t\tH0 [km/s/Mpc]\tw0\t\twa')
	print(f'DESI+CMB+Pantheon:')
	Omega_m = 0.3085
	dOmega_m = 0.0068
	H0 = 68.03
	dH0 = 0.72
	w0 = -0.827
	dw0 = 0.063
	wa = -0.75
	dwap=0.29
	dwam=0.25
	print(f'{Omega_m:5.4f}+/-{dOmega_m:5.4f}\t{H0:3.2f}+/-{dH0:3.2f}\t{w0:4.3f}+/-{dw0:4.3f}\t{wa:5.4f}+{dwap:3.2f}-{dwam:3.2f}')
	print(f'DESI+CMB+Union3:')
	Omega_m = 0.3230
	dOmega_m = 0.0095
	H0 = 66.53
	dH0 = 0.94
	w0 = -0.65
	dw0 = 0.10
	wa = -1.27
	dwap=0.40
	dwam=0.34
	print(f'{Omega_m:5.4f}+/-{dOmega_m:5.4f}\t{H0:3.2f}+/-{dH0:3.2f}\t{w0:4.3f}+/-{dw0:4.3f}\t{wa:5.4f}+{dwap:3.2f}-{dwam:3.2f}')

	print(f'DESI+CMB+DESY5:')
	Omega_m = 0.3160
	dOmega_m = 0.0065
	H0 = 67.24
	dH0 = 0.66
	w0 = -0.727
	dw0 = 0.067
	wa = -1.05
	dwap=0.31
	dwam=0.27
	print(f'{Omega_m:5.4f}+/-{dOmega_m:5.4f}\t{H0:3.2f}+/-{dH0:3.2f}\t{w0:4.3f}+/-{dw0:4.3f}\t{wa:5.4f}+{dwap:3.2f}-{dwam:3.2f}')


def rho(z,args=None):
	if(args is not None):
		H0       = args.H0 
		Omega_m  = args.Omega_m
		Omega_k  = args.Omega_k
		Omega_r  = args.Omega_r
		Omega_de = args.Omega_de
		w0       = args.w0
		wa       = args.wa
	else:
		H0      = 70 #km/s/Mpc
		Omega_m = 0.3
		Omega_k = 0
		Omega_r = 4.203069657615112e-05 / (H0/100)**2
		Omega_de = 1.0 - (Omega_m+Omega_k+Omega_r)



	a = 1./(1.+z)
	o_r  = Omega_r*(1+z)**4
	o_m  = Omega_m*(1+z)**3
	o_k  = Omega_k*(1+z)**2
	o_de = Omega_de*(a**(-3*(1+w0+wa)))*np.exp(-3*wa*(1-a))

	#rho in terms of the critical density
	return o_r+o_m+o_k+o_de

def H(z,args=None):

	if(args is not None):
		H0       = args.H0 
		Omega_m  = args.Omega_m
		Omega_k  = args.Omega_k
		Omega_r  = args.Omega_r
		Omega_de = args.Omega_de
		w0       = args.w0
		wa       = args.wa
		#print('here')
	else:
		H0      = 70 #km/s/Mpc
		Omega_m = 0.3
		Omega_k = 0
		Omega_r = 4.203069657615112e-05 / (H0/100)**2
		Omega_de = 1.0 - (Omega_m+Omega_k+Omega_r)

		#LCDM
		w0       = -1.0
		wa       = 0.0
		#https://arxiv.org/pdf/2404.03002
		#DESI + CMB+Pantheon
		w0       = -0.827
		wa       = -0.75
		#https://arxiv.org/pdf/2404.03002
		#DESI + CMB+Union3
		w0       = -0.64
		wa       = -1.27


		#https://arxiv.org/pdf/2404.03002
		#DESI + CMB+DESY5
		w0       = -0.727
		wa       = -1.05

	a = 1./(1.+z)
	o_r  = Omega_r*(1+z)**4
	o_m  = Omega_m*(1+z)**3
	o_k  = Omega_k*(1+z)**2
	o_de = Omega_de*(a**(-3*(1+w0+wa)))*np.exp(-3*wa*(1-a))

	#return hubble parameter in km/s/Mpc
	return H0*np.sqrt(o_r+o_m+o_k+o_de)

def H_a(a,args=None):
	z = 1./a-1

	#return Hubble parameter in km/s/Mpc
	return H(z)

def dtdz(z, args=None):
	Hz = H(z,args=args)
	return 1./((1+z)*Hz)

def dtdlnz(lnz, args=None):
	z = np.exp(lnz)
	Hz = H(z,args=args)
	return z*1./((1+z)*Hz)

def t(z, args=None):

#	tint, _ = quad(dtdz,0,z,args=args)

	if(isinstance(z,float)):

		#scalar z
		if(z<np.exp(-100)):
			return 0
		lnz = np.log(z)
		tint, _ = quad(dtdlnz,-100,lnz,args=args,limit=100,epsabs=1.0e-10,epsrel=1.0e-10,maxp1=100,limlst=100)
	else:
		#array
		xi = np.where(z<np.exp(-100))[0]
		zz = z.copy()
		zz[xi] = np.exp(-100)
		lnz = np.log(zz)
		tint = np.zeros_like(z)
		for i in range(len(z)):
			tint[i], _ = quad(dtdlnz,-100,lnz[i],args=args,limit=100,epsabs=1.0e-10,epsrel=1.0e-10,maxp1=100,limlst=100)

	#conver Mpc * s / km to Gyr
	mpc_in_km = 3.26156 * 1.0e6 * 31536000 * 2.99792458e5
	gyr_per_s = 1./(1.0e9 * 31536000)
	k_H_inv = mpc_in_km * gyr_per_s 
	#return lookback time
	#at redshift z in Gyr
	return k_H_inv * tint

def H0_in_inv_gyr(args=None):
	if(args is not None):
		H0       = args.H0 
	else:
		H0      = 70 #km/s/Mpc

	print(f'H0 = {H0}')

	#conver Mpc * s / km to Gyr
	mpc_in_km = 3.26156 * 1.0e6 * 31536000 * 2.99792458e5
	gyr_per_s = 1./(1.0e9 * 31536000)
	k_H_inv = mpc_in_km * gyr_per_s 
	return H0/k_H_inv

def get_redshift_array(n,zmax=1000):

    z = np.zeros(n+1)
    z[0] = 0
    z[1:] = np.exp(np.linspace(np.log(0.01),np.log(zmax),n))

    return z

def main():


    #create the command line argument parser
    parser = create_parser()

    #store the command line arguments
    args   = parser.parse_args()

    #check for special known cosmologies
    args = special_cosmologies(args)

    #adjust omega radiation and de 
    args = adjust_omegas(args)

    #print cosmology
    if(args.verbose):
    	#print_reference_histories()
    	print_cosmology(args)

    H0_inv = 1./H0_in_inv_gyr(args=args)
    print(f'1/H0 in gyr = {H0_inv}.')

    z = 0.01
    tl = t(z,args=args)
    print(f'Time to redshift z = {z} is t = {tl} Gyr.')
    z = 1.
    tl = t(z,args=args)
    print(f'Time to redshift z = {z} is t = {tl} Gyr.')
    z = 3.
    tl = t(z,args=args)
    print(f'Time to redshift z = {z} is t = {tl} Gyr.')
    z = 15.
    tl = t(z,args=args)
    print(f'Time to redshift z = {z} is t = {tl} Gyr.')
    z = 1e12
    tl = t(z,args=args)
    print(f'Time to redshift z = {z} is t = {tl} Gyr.')


    z = 2.0
    tl = t(z,args=args)
    print(f'Time to redshift z = {z} is H0t = {tl/H0_inv} Gyr.')
    z = 4.
    tl = t(z,args=args)
    print(f'Time to redshift z = {z} is H0t = {tl/H0_inv} Gyr.')
    z = 6.
    tl = t(z,args=args)
    print(f'Time to redshift z = {z} is H0t = {tl/H0_inv} Gyr.')
    z = 8.
    tl = t(z,args=args)
    print(f'Time to redshift z = {z} is H0t = {tl/H0_inv} Gyr.')
    z = 10.
    tl = t(z,args=args)
    print(f'Time to redshift z = {z} is H0t = {tl/H0_inv} Gyr.')
    z = 12.
    tl = t(z,args=args)
    print(f'Time to redshift z = {z} is H0t = {tl/H0_inv} Gyr.')
    z = 14.
    tl = t(z,args=args)
    print(f'Time to redshift z = {z} is H0t = {tl/H0_inv} Gyr.')
    z = 16.
    tl = t(z,args=args)
    print(f'Time to redshift z = {z} is H0t = {tl/H0_inv} Gyr.')

    tuni = t(1e12,args=args)
    z = 2.0
    tl = t(z,args=args)
    print(f'Time to redshift z = {z} is t = {tuni-tl} Gyr.')
    z = 2.0
    tl = t(z,args=args)
    print(f'Time to redshift z = {z} is t = {tuni-tl} Gyr.')
    z = 4.
    tl = t(z,args=args)
    print(f'Time to redshift z = {z} is t = {tuni-tl} Gyr.')
    z = 6.
    tl = t(z,args=args)
    print(f'Time to redshift z = {z} is t = {tuni-tl} Gyr.')
    z = 8.
    tl = t(z,args=args)
    print(f'Time to redshift z = {z} is t = {tuni-tl} Gyr.')
    z = 10.
    tl = t(z,args=args)
    print(f'Time to redshift z = {z} is t = {tuni-tl} Gyr.')
    z = 12.
    tl = t(z,args=args)
    print(f'Time to redshift z = {z} is t = {tuni-tl} Gyr.')
    z = 14.
    tl = t(z,args=args)
    print(f'Time to redshift z = {z} is t = {tuni-tl} Gyr.')
    z = 16.
    tl = t(z,args=args)
    print(f'Time to redshift z = {z} is t = {tuni-tl} Gyr.')


    z = get_redshift_array(1000)
    a = 1./(1+z)
    tl = t(z,args=args)
    tlb = tl/H0_inv

    np.savetxt("a_vs_H0t.txt",np.asarray([tlb,a]).T)
    plot_scale_factor_vs_look_back_time(tlb,a,args=args)

if __name__=="__main__":
	main()
