import numpy as np
import h5py
import sys
import os
import copy
import os.path
import pickle as pickle
import pkg_resources
import multiprocessing
import types
from multiprocessing import Process, Queue
from time import sleep

from scipy import interpolate
from scipy.io import FortranFile
import scipy.integrate as integrate


from src.fortran.get_ptcl_py import get_ptcl_py
from src.fortran.get_cell_py import get_cell_py
from src.fortran.get_amr_py import get_amr_py
from src.fortran.find_domain_py import find_domain_py


class veluga:

	def __init__(self, header, num_thread=1):

		##-----
		## INIT Setting
		##-----
		self.sun_met 	= 0.02
		self.num_thread = int(num_thread)
		self.r_gal_tname = None

		void 	= self.rdheader(header)

##-----
## SIMPLE FTNS
##-----

	def rdheader(self, fname):
		fdata = open(fname, 'r')
		lines = fdata.readlines()

		header 	= types.SimpleNamespace()
		for line in lines:
			
			if '#' in line: continue
			if 'pp_' in line: continue
			if 'makebr_' in line: continue
			
			s0 	= line.split(';')[0]
			exec(s0, header.__dict__)

		fdata.close()

		self.header = header

	def allocate(self, npart, type=None):

		if(type=='part'):
			dtype = [('xx', '<f8'), ('yy', '<f8'), ('zz', '<f8'),
				('vx', '<f8'), ('vy', '<f8'), ('vz', '<f8'),
				('mp', '<f8'), ('ap', '<f8'), ('zp', '<f8'), ('gyr', '<f8'), ('sfact', '<f8'), ('redsh', '<f8'), ('family', np.int32) , ('domain', np.int32), 
				('id', np.int64), ('KE', '<f8'), ('UE', '<f8'), ('PE', '<f8')]
		elif(type=='cell'):
			print('add hydrovariables here')
			dtype = [('xx', '<f8'), ('yy', '<f8'), ('zz', '<f8'),
				('vx', '<f8'), ('vy', '<f8'), ('vz', '<f8'),
				('level', np.int32), ('dx', '<f8'), ('den', '<f8'), ('temp', '<f8'), ('zp', '<f8'), ('mp', '<f8'), 
				('KE', '<f8'), ('UE', '<f8'), ('PE', '<f8')]

		return np.zeros(npart, dtype=dtype)

	def errorout(self, str):
		print('%123123----- VELUGA -----')
		print(' ')
		for str0 in [str]: print('		' + str)
		print(' ')
		print('%123123----------')

	def stat(self, array):
		print('MIN = ', np.amin(array))
		print('MAX = ', np.amax(array))
		print('AVG = ', np.average(array))
		print('MED = ', np.median(array))

	def ramses_part_input(self, ramobj):

		ramobj.n_thread = np.int32(self.num_thread)

		if(self.header.idtype == 'long'): ramobj.r_idtype = np.int32(-1)
		elif(self.header.idtype == 'long64'): ramobj.r_idtype = np.int32(1)
		else: self.errout('Wrong value for idtype')

		if(self.header.famtype == 'old'): ramobj.r_famtype = np.int32(-1)
		elif(self.header.famtype == 'new'): ramobj.r_famtype = np.int32(1)
		else: self.errout('Wrong value for famtype')

		ramobj.r_skip_domain = np.int32(-1)
		if(self.header.skiprd_domain == 1): ramobj.r_skip_domain = np.int32(1)

		ramobj.r_skip_time = np.int32(-1)
		if(self.header.skiprd_time == 1): ramobj.r_skip_time = np.int32(1)

		ramobj.r_skip_metal = np.int32(-1)
		if(self.header.skiprd_metal == 1): ramobj.r_skip_metal = np.int32(1)

		ramobj.dir_raw = self.header.dir_raw.ljust(1000)

	def ramses_cell_input(self, ramobj):

		info = self.g_info(1)

		ramobj.n_thread = np.int32(self.num_thread)
		ramobj.n_mpi = np.int32(self.header.ndomain)
		ramobj.n_dim = np.int32(3)
		ramobj.levmax = np.int32(info['levmax'])
		ramobj.dir_raw = self.header.dir_raw.ljust(1000)

	def ramses_amr_input(self, ramobj):

		info = self.g_info(1)

		ramobj.n_thread = np.int32(self.num_thread)
		ramobj.n_mpi = np.int32(self.header.ndomain)
		ramobj.n_dim = np.int32(3)
		ramobj.levmax = np.int32(info['levmax'])
		ramobj.dir_raw = self.header.dir_raw.ljust(1000)

##-----
## GET FTNS
##-----
	##-----
	## Get Ramses Info
	##-----
	def g_info(self, snapnum):
	
		fname   = self.header.dir_raw + "output_%0.5d"%snapnum + "/info_%0.5d"%snapnum + ".txt"
	   
		fdata1  = np.loadtxt(fname, dtype=object, max_rows=6, delimiter='=')
		fdata2  = np.loadtxt(fname, dtype=object, skiprows=7, max_rows=11, delimiter='=')

		kpc     = np.double(3.086e21) #3.08568025e21
		twopi   = np.double(6.2831853e0)
		hplanck = np.double(6.6262000e-27)
		eV      = np.double(1.6022000e-12)
		kB      = np.double(1.3806200e-16)
		clight  = np.double(2.9979250e+10)
		Gyr     = np.double(3.1536000e+16)
		X       = np.double(0.76)
		Y       = np.double(0.24 )
		rhoc    = np.double(1.8800000e-29)
		mH      = np.double(1.6600000e-24)
		mu_mol  = np.double(1.2195e0)
		G       = np.double(6.67259e-8)
		m_sun   = np.double(1.98892e33)

		cgs 	= {'kpc':kpc, 'hplanck':hplanck, 'eV':eV, 'kB':kB, 'clight':clight, 'Gyr':Gyr, 'mH':mH, 'G':G, 'm_sun':m_sun}

		info    = {'ncpu':0, 'ndim':0, 'levmin':0, 'levmax':0, 'aexp':0, 
				'H0':0, 'oM':0, 'oB':0, 'oL':0, 'unit_l':0, 'unit_d':0, 'unit_T2':0, 'nH':0, 
				'unit_t':0, 'kms':0, 'unit_m':0, 'hindex':0, 'cgs':cgs}

		info['ncpu']    = np.int32(fdata1[0][1])
		info['ndim']    = np.int32(fdata1[1][1])
		info['levmin']  = np.int32(fdata1[2][1])
		info['levmax']  = np.int32(fdata1[3][1])
		info['aexp']    = np.double(fdata2[2][1])
		info['H0']      = np.double(fdata2[3][1])
		info['oM']      = np.double(fdata2[4][1])
		info['oL']      = np.double(fdata2[5][1])
		info['oB']      = np.double(fdata2[7][1])
		info['unit_l']  = np.double(fdata2[8][1])
		info['unit_d']  = np.double(fdata2[9][1])
		info['unit_t']  = np.double(fdata2[10][1])
		info['unit_T2'] = np.double(1.66e-24) / np.double(1.3806200e-16) * np.double(info['unit_l'] / info['unit_t'])**2
		info['nH']      = np.double(0.76) / np.double(1.66e-24) * info['unit_d']
		info['kms']     = info['unit_l'] / info['unit_t'] / 1e5
		info['unit_m']  = info['unit_d'] * info['unit_l']**3

		info['hindex']  = np.double(np.loadtxt(fname, dtype=object, skiprows=21)[:,1:])

		return info		

	##-----
	## Get Physical time from conformal time
	##-----
	def g_ptime(self, n_snap, t_conf):

		#----- Initial settings
		info 	= self.g_info(n_snap)
		

		#----- Allocate
		time 	= np.zeros(len(t_conf), dtype=[('sfact', '<f8'), ('gyr', '<f8'), ('redsh', '<f8')])

		#----- Get Conformal Time - Scale factor table
		c_table 	= self.g_cfttable(info)

		#----- Interpolation
		lint = interpolate.interp1d(c_table['conft'],c_table['sfact'],kind = 'quadratic')
		time['sfact'][:] 	= lint(t_conf)
		time['redsh'][:]	= 1./time['sfact'][:] - 1.

		#----- Get Gyr from scale factor table
		g_table = self.g_gyrtable(info)

		#----- Interplation
		lint = interpolate.interp1d(g_table['redsh'],g_table['gyr'],kind = 'quadratic')
		t0  = lint( 1./info['aexp'] - 1.)

		time['gyr'][:] 		= lint(time['redsh'][:]) - t0

		return time

	##-----
	## Generate or Load Confal-Gyr Table
	##-----
	def g_cfttable_ftn(self, X, oM, oL):
		return 1./(X**3 * np.sqrt(oM/X**3 + oL))

	def g_cfttable(self, info):


		H0     = info['H0']
		oM     = info['oM']
		oL     = info['oL']

		# reference path is correct?
		fname 	= 'table/cft_%0.5d'%(H0*1000.) + '_%0.5d'%(oM*100000.) + '_%0.5d'%(oL*100000.) + '.pkl'
		isfile 	= os.path.isfile(fname)

		if(isfile==True):
			with open(fname, 'rb') as f:
				data = pickle.load(f)
		else:
			n_table = np.int32(10000)
			data    = np.zeros(n_table, dtype=[('sfact','<f8'), ('conft','<f8')])
			data['sfact'][:]    = np.array(range(n_table),dtype='<f8')/(n_table - 1.) * 0.98 + 0.02

			ind     = np.array(range(n_table),dtype='int32')
			for i in ind:
				data['conft'][i]  = integrate.quad(self.g_cfttable_ftn,data['sfact'][i],1.,args=(oM,oL))[0] * (-1.)

			with open(fname, 'wb') as f:
				pickle.dump(data, f, pickle.HIGHEST_PROTOCOL)

		return data

	##-----
	## Generate or Load Sfactor-Gyr Table
	##-----
	def g_gyrtable_ftn(self, X, oM, oL):
		return 1./(1.+X)/np.sqrt(oM*(1.+X)**3 + oL)

	def g_gyrtable(self, info):

		H0     = info['H0']
		oM     = info['oM']
		oL     = info['oL']

		fname   = 'table/gyr_%0.5d'%(H0*1000.) + '_%0.5d'%(oM*100000.) + '_%0.5d'%(oL*100000.) + '.pkl'
		isfile = os.path.isfile(fname)

		if(isfile==True):
			with open(fname, 'rb') as f:
				data = pickle.load(f)
		else:
			n_table = np.int32(10000)
			data    = np.zeros(n_table, dtype=[('redsh','<f8'),('gyr','<f8')])
			data['redsh'][:]    = 1./(np.array(range(n_table),dtype='<f8')/(n_table - 1.) * 0.98 + 0.02) - 1.
			data['gyr'][0]  = 0.

			ind     = np.array(range(n_table),dtype='int32')
			for i in ind:
				data['gyr'][i]  = integrate.quad(self.g_gyrtable_ftn,0.,data['redsh'][i], args=(oM, oL))[0]
				data['gyr'][i]  *= (1./H0 * np.double(3.08568025e19) / np.double(3.1536000e16))

			with open(fname, 'wb') as f:
				pickle.dump(data, f, pickle.HIGHEST_PROTOCOL)

		return data


	##-----
	## Get Domain containing a volume dl^3 centered at (x0, y0, z0)
	##	x0, y0, z0, dl should be given in pKpc unit
	##-----
	def g_domain(self, n_snap, x0, y0, z0, dl):
		info 	= self.g_info(n_snap)

		x0 = np.array(x0, dtype='<f8')
		y0 = np.array(y0, dtype='<f8')
		z0 = np.array(z0, dtype='<f8')
		dl = np.array(dl, dtype='<f8')

		n_gal 	= np.size(x0)
		n_mpi 	= info['ncpu']

		if(np.size(x0) != np.size(dl)): dl = x0*0. + np.amax(dl)

		x0_s 	= x0 * (info['cgs']['kpc'] / info['unit_l'])
		y0_s 	= y0 * (info['cgs']['kpc'] / info['unit_l'])
		z0_s 	= z0 * (info['cgs']['kpc'] / info['unit_l'])
		dl_s 	= dl * (info['cgs']['kpc'] / info['unit_l'])

		#----- Fortran Settings
		find_domain_py.n_ptcl 	= np.int32(n_gal)
		find_domain_py.n_mpi 	= np.int32(n_mpi)
		find_domain_py.n_thread	= np.int32(self.num_thread)
		find_domain_py.n_dim 	= np.int32(3)
		find_domain_py.levmax	= np.int32(info['levmax'])


		##----- Find Domain
		find_domain_py.find_domain(x0_s, y0_s, z0_s, dl_s, info['hindex'])
		dom_list 	= find_domain_py.dom_list

		if(np.size(np.where(dom_list > 0)) == 0):
			self.errorout('No domains are found (error)')
			find_domain_py.find_domain_free()
			return -1


		##----- Merge
		dom_all 	= np.zeros(n_gal*n_mpi, dtype=np.int32) - 1

		i0 = 0
		for i in range(0, n_gal):
			cut = np.array(np.where(dom_list[i,:] > 0),dtype=np.int32) + 1

			if(np.size(cut) == 0):
				print(' ? ')
				return -1

			dom_all[i0:i0+np.size(cut)] = cut
			i0 = i0 + np.size(cut)

		dom_all = np.unique(dom_all)
		return dom_all[1:]

	##-----
	## Get Particle within a volume dl^3 centered at (x0, y0, z0)
	##  x0, y0, z0, dl should be given in pKpc unit
	##  ptype as 'all' 'dm' or 'star'
	##
	##	if dom_list is given as an numpy integer array, it overwries the domain found with the given volume and return all ptcls in the argued domain
	##-----
	def g_part(self, n_snap, x0, y0, z0, dl, ptype='all', dom_list=None, g_simunit=False, g_ptime=False):

		##----- Initial Settings

		info = self.g_info(n_snap)
		dmp_mass    = 1.0/(self.header.neff*self.header.neff*self.header.neff)*(info['oM'] - info['oB'])/info['oM']

		if(dom_list is None):
			dom_list = self.g_domain(n_snap, x0, y0, z0, dl)
			x0_s = np.double(x0 * (info['cgs']['kpc'] / info['unit_l']))
			y0_s = np.double(y0 * (info['cgs']['kpc'] / info['unit_l']))
			z0_s = np.double(z0 * (info['cgs']['kpc'] / info['unit_l']))
			dl_s = np.double(dl * (info['cgs']['kpc'] / info['unit_l']))
		else:
			x0_s = np.double(0.5)
			y0_s = np.double(0.5)
			z0_s = np.double(0.5)
			z0_s = np.double(0.5)
			dl_s = np.double(1e8)

		##----- Fortran Settigns
		self.ramses_part_input(get_ptcl_py)
		if(ptype == 'all'): get_ptcl_py.r_ptype = np.int32(0)
		elif(ptype == 'dm'): get_ptcl_py.r_ptype = np.int32(1)
		elif(ptype == 'star'): get_ptcl_py.r_ptype = np.int32(2)
		get_ptcl_py.dmp_mass = dmp_mass

		## Get by Fortran
		get_ptcl_py.get_ptcl_box(n_snap, x0_s, y0_s, z0_s, dl_s, dom_list)

		## Allocate
		n_ptcl = np.size(get_ptcl_py.p_dbl[:,0])
		part 	= self.allocate(n_ptcl, type='part')

		part['xx'] 	= get_ptcl_py.p_dbl[:,0]
		part['yy'] 	= get_ptcl_py.p_dbl[:,1]
		part['zz'] 	= get_ptcl_py.p_dbl[:,2]
		part['vx'] 	= get_ptcl_py.p_dbl[:,3]
		part['vy'] 	= get_ptcl_py.p_dbl[:,4]
		part['vz'] 	= get_ptcl_py.p_dbl[:,5]
		part['mp'] 	= get_ptcl_py.p_dbl[:,6]
		part['ap'] 	= get_ptcl_py.p_dbl[:,7]
		part['zp'] 	= get_ptcl_py.p_dbl[:,8]
		part['id'] 	= get_ptcl_py.p_lnt[:,0]
		part['family'] = get_ptcl_py.p_lnt[:,1]

		get_ptcl_py.get_ptcl_deallocate()

		##
		if(g_simunit==False):
			part['xx'] 	*= (info['unit_l'] / info['cgs']['kpc'])
			part['yy'] 	*= (info['unit_l'] / info['cgs']['kpc'])
			part['zz'] 	*= (info['unit_l'] / info['cgs']['kpc'])
			part['vx'] 	*= (info['kms'])
			part['vy'] 	*= (info['kms'])
			part['vz'] 	*= (info['kms'])
			part['mp']	*= (info['unit_m'] / info['cgs']['m_sun'])

		if(g_ptime==True):
			ptime = self.g_ptime(n_snap, part['ap'])
			part['gyr'] = ptime['gyr']
			part['sfact'] = ptime['sfact']
			part['redsh'] = ptime['redsh']


		return part

	##-----
	## Get Cell within a volume dl^3 centered at (x0, y0, z0)
	##  x0, y0, z0, dl should be given in pKpc unit
	##
	##
	##	if dom_list is given as an numpy integer array, it overwries the domain found with the given volume and return all cells in the argued domain
	##-----
	def g_cell(self, n_snap, x0, y0, z0, dl, dom_list=None, g_simunit=False):

		##----- Initial Settings

		info = self.g_info(n_snap)
		dmp_mass    = 1.0/(self.header.neff*self.header.neff*self.header.neff)*(info['oM'] - info['oB'])/info['oM']

		if(dom_list is None):
			dom_list = self.g_domain(n_snap, x0, y0, z0, dl)
			x0_s = np.double(x0 * (info['cgs']['kpc'] / info['unit_l']))
			y0_s = np.double(y0 * (info['cgs']['kpc'] / info['unit_l']))
			z0_s = np.double(z0 * (info['cgs']['kpc'] / info['unit_l']))
			dl_s = np.double(dl * (info['cgs']['kpc'] / info['unit_l']))
		else:
			dom_list = np.int32(dom_list)
			x0_s = np.double(0.5)
			y0_s = np.double(0.5)
			z0_s = np.double(0.5)
			z0_s = np.double(0.5)
			dl_s = np.double(1e8)

		##----- Fortran Settigns
		self.ramses_cell_input(get_cell_py)

		

		##----- Get cell by Fortran
		get_cell_py.get_cell_box(n_snap, x0_s, y0_s, z0_s, dl_s, dom_list)


		##----- Allocate
		n_cell = np.size(get_cell_py.p_dbl[:,0])
		cell 	= self.allocate(n_cell, type='cell')

		##----- Get Hydro descriptor
		cell['xx'] 	= get_cell_py.p_dbl[:,0]
		cell['yy'] 	= get_cell_py.p_dbl[:,1]
		cell['zz'] 	= get_cell_py.p_dbl[:,2]

		hind = np.int32(3)
		for htag in self.header.hydro_variables:
			if htag == 'skip': continue
			cell[htag] = get_cell_py.p_dbl[:,hind]
			hind += 1

		cell['dx'] 	= get_cell_py.p_dbl[:,hind+1]
		cell['level'] = get_cell_py.p_int[:,0]

		
		##----- Memory free
		get_cell_py.get_cell_deallocate()

		##
		cell['mp'] = cell['den'] * info['unit_d'] * (cell['dx']*info['unit_l'])**3 / info['cgs']['m_sun'] ## Msun
		cell['UE'] = cell['temp']/cell['den']/(5./3.-1.) * info['unit_T2'] / (1.66e-24) * 1.38049e-23 * 1e-3 ## [km/s]^2
		cell['temp']*= (1./cell['den'])

		if(g_simunit==False):
			cell['xx'] 	*= (info['unit_l'] / info['cgs']['kpc'])
			cell['yy'] 	*= (info['unit_l'] / info['cgs']['kpc'])
			cell['zz'] 	*= (info['unit_l'] / info['cgs']['kpc'])
			cell['vx'] 	*= (info['kms'])
			cell['vy'] 	*= (info['kms'])
			cell['vz'] 	*= (info['kms'])
			cell['dx']	*= (info['unit_l'] / info['cgs']['kpc'])
			cell['den'] *= info['nH']		# [/cc]
			cell['temp']*= info['unit_T2'] 		# [K/mu]

		

		return cell

	##-----
	## Get AMR within a volume dl^3 centered at (x0, y0, z0)
	##	output includes x, y, z, dx (not hydro variables)
	##	active for non-hydro runs
	##  x0, y0, z0, dl should be given in pKpc unit
	##
	##
	##	if dom_list is given as an numpy integer array, it overwries the domain found with the given volume and return all cells in the argued domain
	##-----
	def g_amr(self, n_snap, x0, y0, z0, dl, dom_list=None, g_simunit=False):

		##----- Initial Settings

		info = self.g_info(n_snap)
		dmp_mass    = 1.0/(self.header.neff*self.header.neff*self.header.neff)*(info['oM'] - info['oB'])/info['oM']

		if(dom_list is None):
			dom_list = self.g_domain(n_snap, x0, y0, z0, dl)
			x0_s = np.double(x0 * (info['cgs']['kpc'] / info['unit_l']))
			y0_s = np.double(y0 * (info['cgs']['kpc'] / info['unit_l']))
			z0_s = np.double(z0 * (info['cgs']['kpc'] / info['unit_l']))
			dl_s = np.double(dl * (info['cgs']['kpc'] / info['unit_l']))
		else:
			dom_list = np.int32(dom_list)
			x0_s = np.double(0.5)
			y0_s = np.double(0.5)
			z0_s = np.double(0.5)
			z0_s = np.double(0.5)
			dl_s = np.double(1e8)

		##----- Fortran Settigns
		self.ramses_amr_input(get_amr_py)

		

		##----- Get cell by Fortran
		get_amr_py.get_amr_box(n_snap, x0_s, y0_s, z0_s, dl_s, dom_list)


		##----- Allocate
		n_cell = np.size(get_amr_py.p_dbl[:,0])
		cell 	= self.allocate(n_cell, type='cell')

		##----- Get Hydro descriptor
		cell['xx'] 	= get_amr_py.p_dbl[:,0]
		cell['yy'] 	= get_amr_py.p_dbl[:,1]
		cell['zz'] 	= get_amr_py.p_dbl[:,2]
		cell['dx'] 	= get_amr_py.p_dbl[:,3]
		cell['level'] = get_amr_py.p_int[:,0]

		##----- Memory free
		get_amr_py.get_amr_deallocate()

		
		if(g_simunit==False):
			cell['xx'] 	*= (info['unit_l'] / info['cgs']['kpc'])
			cell['yy'] 	*= (info['unit_l'] / info['cgs']['kpc'])
			cell['zz'] 	*= (info['unit_l'] / info['cgs']['kpc'])
			cell['dx']	*= (info['unit_l'] / info['cgs']['kpc'])

		return cell
##-----
## LOAD CATALOG ROUTINES
##-----

	##-----
	## Load Galaxy
	##-----
	class r_gal_parallel:
		def __init__(self, vrobj, h5data, gidlist, galdata):
			self.galdata = galdata
			self.gidlist = gidlist
			self.h5data = h5data
			self.vrobj = vrobj

		def r_gal_input(self, start, end):
			for i in range(start, end):
				idstr  = 'ID_%0.6d'%self.gidlist[i]

				for name in self.vrobj.r_gal_tname:

					if(name=='snapnum' or name=='redsh' or name=='Aexp'): continue
					if(name=='isclump' or name=='Domain_List'): str0 = idstr + '/' + name
					else: str0 = idstr + '/G_Prop/G_' + name
					xdata = self.h5data.get(str0)
					self.galdata[name][i] = np.array(xdata)

				#for name in self.vrobj.gal_prop:
				#	xdata = self.h5data.get(idstr + "/G_Prop/G_" + name)
				#	self.galdata[name][i] = np.array(xdata)

				#for name in self.vrobj.vr_general:
				#	xdata = self.h5data.get("/" + name)
				#	self.galdata[name][i] = np.array(xdata)
				#for name in self.vrobj.vr_galinfo:
				#	if(name!='snapnum'): 
				#		xdata = self.h5data.get(idstr + "/" + name)
				#		self.galdata[name][i] = np.array(xdata)


		def r_gal_input_p(self, start, end, q):
			for i in range(start, end):
				idstr  = 'ID_%0.6d'%self.gidlist[i]

				for name in self.vrobj.r_gal_tname:

					if(name=='snapnum' or name=='redsh' or name=='Aexp'): continue
					if(name=='isclump' or name=='Domain_List'): str0 = idstr + '/' + name
					else: str0 = idstr + '/G_Prop/G_' + name

					xdata = self.h5data.get(str0)
					self.galdata[name][i] = np.array(xdata)

				#for name in self.vrobj.vr_galprop:
				#	xdata = self.h5data.get(idstr + "/G_Prop/G_" + name)
				#	self.galdata[name][i] = np.array(xdata)
	
				#for name in self.vrobj.vr_general:
				#	xdata = self.h5data.get("/" + name)
				#	self.galdata[name][i] = np.array(xdata)
	
				#for name in self.vrobj.vr_galinfo:
				#	if(name!='snapnum'):
				#		xdata = self.h5data.get(idstr + "/" + name)
				#		self.galdata[name][i] = np.array(xdata)

			q.put((start, end, self.galdata[start:end]))

	def r_gal(self, n_snap, id0, mrange=None, horg='g'):
		"""
		Load Galaxy/Halo Catalog Data.

		This method retrieves galaxy or halo catalog data for a given snapshot and object ID.
		
		Parameters
		----------
		n_snap : int
			Snapshot number
		id0 : int
			Object ID. Use a negative value to retrieve all objects in the snapshot.
		horg : {'g', 'h'}
			A flag to specify the object type. Use 'g' for galaxies and 'h' for halos.
			Default is 'g'.

		Returns
		-------
		structured_array
			A structured array containing information about the objects.
		
		Raises
		------
		ValueError
			If `n_snap` is not a valid snapshot number.
			If `horg` is not one of {'g', 'h'}.
	
		Examples
		--------
		>>> g = r_gal(100, -1, horg='g')
		>>> print(g['id'])
		>>> h = r_gal(200, 1, mrange=(1e10, 1e12), horg='h')
		>>> print(h['mass'])

		"""

		# Path setting
		if(horg=='h'): fname = self.header.dir_catalog + 'Halo/VR_Halo/snap_%0.4d'%n_snap + '.hdf5'
		elif(horg=='g'): fname = self.header.dir_catalog + 'Galaxy/VR_Galaxy/snap_%0.4d'%n_snap+'.hdf5'
		else:
			self.errout([' Wrong argument for the horg ', '     horg = "g" (for galaxy) or "h" (for halo)'])


		# Open hdf5
		dat     = h5py.File(fname, 'r')

		# Get Some bulk first
		flux_list 	= np.array(dat.get('Flux_List'), dtype='str')
		sfr_r 		= np.array(dat.get('SFR_R'), dtype='<f8')
		sfr_t 		= np.array(dat.get('SFR_T'), dtype='<f8')
		mag_r 		= np.array(dat.get('MAG_R'), dtype='<f8')
		conf_r 		= np.array(dat.get('CONF_R'), dtype='<f8')
		mass_tot 	= np.array(dat.get('Mass_tot'), dtype='<f8')
		r_halfmass	= np.array(dat.get('R_HalfMass'), dtype='<f8')
		mvir 		= np.array(dat.get('Mvir'), dtype='<f8')
		rvir 		= np.array(dat.get('Rvir'), dtype='<f8')
		idlist 		= np.array(dat.get('ID'), dtype=np.int32)
		conf_m 		= np.array(dat.get('CONF_M'), dtype='<f8')
		conf_n 		= np.array(dat.get('CONF_N'), dtype='<f8')
		aexp 		= np.array(dat.get('ID_000001/Aexp'), dtype='<f8')

		# Set column list
		dtype=[]

		##----- Original catalog
		for name in self.header.column_list:
			if(name=='ID' or name=='hostHaloID' or name=='numSubStruct' or name=='Structuretype' or name=='npart'): dtype=dtype+[(name, np.int32, (1,))]
			elif(name=='ID_mbp'): dtype=dtype+[(name, np.int64, (1,))]
			else: dtype=dtype+[(name, '<f8', (1,))]

		##----- Gal prop with post-processing
		for name in self.header.gal_prop:
			name0 = name
			if(name0 == 'sfr'): name0 = 'SFR'
			elif(name0 == 'abmag'): name0 = 'ABmag'
			elif(name0 == 'sb'): name0 = 'SB'
			elif(name0 == 'conf' or name0 == 'confrac'): name0 = 'CONF'

			if(name0=='SFR'):
				dtype = dtype + [(name0, '<f8', (np.size(sfr_r),))]
			if(name0=='ABmag' or name0 == 'SB'):
				for flist in flux_list:
					name1 = name0 + '_' + flist
					dtype = dtype + [(name1, '<f8', (np.size(mag_r),))]

			if(name0=='CONF'):
				dtype = dtype + [('ConFrac_N','<f8', (np.size(conf_r),)), ('ConFrac_M','<f8', (np.size(conf_r),))]


		##----- Other properties
		dtype = dtype + [('snapnum','<i4', (1,)), ('redsh','<f8', (1,)), ('Aexp', '<f8', (1,)), ('Domain_List', np.int32, (self.header.ndomain,)), ('isclump','<i4', (1,))]
		dtype2 = [('flux_list', 'U10', (np.size(flux_list),)), ('CONF_R','<f8', (np.size(conf_r),)), ('MAG_R', '<f8', (np.size(mag_r),)), ('SFR_R','<f8', (np.size(sfr_r),)), ('SFR_T','<f8', (np.size(sfr_t),))]

		
		#if(self.r_gal_tname == None): 
		self.r_gal_tname = [item[0] for item in dtype]

		# Extract
		if id0 > 0:
			idlist = idlist[idlist == id0]
		else:
			if not (mrange is None):
				if horg == 'g': idlist = idlist[(mass_tot >= mrange[0]) * (mass_tot < mrange[1])]
				if horg == 'h': idlist = idlist[(mvir >= mrange[0]) * (mvir < mrange[1])]

		n_gal   = len(idlist)

		##----- ALLOCATE
		galdata=np.zeros(n_gal, dtype=dtype+dtype2)

		galdata['snapnum'] = n_snap
		galdata['Aexp'] = aexp
		galdata['redsh'] = 1./aexp - 1.

		galdata['flux_list'] = flux_list
		galdata['CONF_R'] = conf_r
		galdata['SFR_R'] = sfr_r
		galdata['SFR_T'] = sfr_t
		galdata['MAG_R'] = mag_r

		#with multiprocessing.Pool(self.num_thread) as pool:
		#    #tasks = [(self.r_gal_input, (i, gid, galdata, dat)) for i, gid in enumerate(idlist)]
		#    void = pool.starmap(self.r_gal_input, [(i, gid, galdata, dat) for i, gid in enumerate(idlist)])
		#    pool.close()
		#    pool.join()
		#if __name__ == '__main__':
		p_input = self.r_gal_parallel(self, dat, idlist, galdata)

		if n_gal < self.num_thread:
			#serially
			p_input.r_gal_input(0, n_gal)
			return p_input.galdata
		else:
			#parallelly (motivated by uhmi.py)
			dind  = np.int32(n_gal / self.num_thread)
			ps = []
			q = Queue()

			for th in range(self.num_thread):
				i0 = th * dind
				i1 = (th+1) * dind
				if(th==0): i0 = 0
				if(th==self.num_thread-1): i1 = n_gal
				p = Process(target=p_input.r_gal_input_p, args=(i0, i1, q))
				ps.append(p)

				p.start()
				while not q.empty():
					i0, i1, dumdata = q.get()
					galdata[i0:i1] = dumdata
			ok = False
			while not ok:
				ok = True
				for idx in np.arange(len(ps)):
					if (ps[idx].is_alive()):
						ok = False
				if(not q.empty()):
					i0, i1, dumdata = q.get()
					galdata[i0:i1] = dumdata
				else:
					sleep(0.5)

			return galdata

	##-----
	## Return cell in a box with dx^3 centered at the galaxy center
	##	
	##-----
	def r_cell(self, n_snap, id0, dx, horg='g'):

		g = self.r_gal(n_snap, id0, horg=horg)

		return g_cell(n_snap, g['Xc'], g['Yc'], g['Zc'], dx)


	##-----
	## Return Member particle ID
	##-----
	def r_pid(self, n_snap, id0, horg='g'):

		# Path setting
		if(horg=='h'): fname = self.header.dir_catalog + 'Halo/VR_Halo/snap_%0.4d'%n_snap + '.hdf5'
		elif(horg=='g'): fname = self.header.dir_catalog + 'Galaxy/VR_Galaxy/snap_%0.4d'%n_snap+'.hdf5'
		else:
			self.errout([' Wrong argument for the horg ', '     horg = "g" (for galaxy) or "h" (for halo)'])

		# Open hdf5
		dat     = h5py.File(fname, 'r')

		str0  = 'ID_%0.6d'%id0 + '/P_Prop/P_ID'
		return np.array(dat.get(str0))

	##-----
	## Return Domain List
	##-----
	def r_domain(self, n_snap, id0, horg='g'):
		# Path setting
		if(horg=='h'): fname = self.header.dir_catalog + 'Halo/VR_Halo/snap_%0.4d'%n_snap + '.hdf5'
		elif(horg=='g'): fname = self.header.dir_catalog + 'Galaxy/VR_Galaxy/snap_%0.4d'%n_snap+'.hdf5'
		else:
			self.errout([' Wrong argument for the horg ', '     horg = "g" (for galaxy) or "h" (for halo)'])

		# Open hdf5
		dat     = h5py.File(fname, 'r')

		str0  = 'ID_%0.6d'%id0 + '/Domain_List'
		dlist0= np.array(dat.get(str0))

		return np.int32(np.where(dlist0 > 0)[0] + 1)

	##-----
	## Load Particle of a galaxy
	##  To do list
	##      *) Halo member load is not implemented
	##
	##	simout=True results in data unit in simulation units
	##-----
	def r_part(self, n_snap, id0, horg='g', g_simunit=False, g_ptime=False):


		## Load requirements
		pid 	= self.r_pid(n_snap, id0, horg=horg)
		dom_list= self.r_domain(n_snap, id0, horg=horg)
		info 	= self.g_info(n_snap)	

		## Initial set
		dmp_mass    = 1.0/(self.header.neff*self.header.neff*self.header.neff)*(info['oM'] - info['oB'])/info['oM']
		n_ptcl 	= np.size(pid)

		## Fortran settings for READ
		self.ramses_part_input(get_ptcl_py)
		if(horg == 'g'): get_ptcl_py.r_ptype = np.int32(2)
		elif(horg == 'h'): get_ptcl_py.r_ptype = np.int32(1)

		if(horg == 'h'): get_ptcl_py.dmp_mass = dmp_mass

		## Get by Fortran
		get_ptcl_py.get_ptcl(n_snap, pid, dom_list)

		## Allocate
		part 	= self.allocate(n_ptcl, type='part')

		part['xx'] 	= get_ptcl_py.p_dbl[:,0]
		part['yy'] 	= get_ptcl_py.p_dbl[:,1]
		part['zz'] 	= get_ptcl_py.p_dbl[:,2]
		part['vx'] 	= get_ptcl_py.p_dbl[:,3]
		part['vy'] 	= get_ptcl_py.p_dbl[:,4]
		part['vz'] 	= get_ptcl_py.p_dbl[:,5]
		part['mp'] 	= get_ptcl_py.p_dbl[:,6]
		part['ap'] 	= get_ptcl_py.p_dbl[:,7]
		part['zp'] 	= get_ptcl_py.p_dbl[:,8]
		part['id'] 	= get_ptcl_py.p_lnt[:,0]
		if(self.header.famtype == 'famtype'): part['family'] = get_ptcl_py.p_lnt[:,1]

		get_ptcl_py.get_ptcl_deallocate()

		##
		if(g_simunit==False):
			part['xx'] 	*= (info['unit_l'] / info['cgs']['kpc'])
			part['yy'] 	*= (info['unit_l'] / info['cgs']['kpc'])
			part['zz'] 	*= (info['unit_l'] / info['cgs']['kpc'])
			part['vx'] 	*= (info['kms'])
			part['vy'] 	*= (info['kms'])
			part['vz'] 	*= (info['kms'])
			part['mp']	*= (info['unit_m'] / info['cgs']['m_sun'])

		if(g_ptime==True):
			ptime = self.g_ptime(n_snap, part['ap'])
			part['gyr'] = ptime['gyr']
			part['sfact'] = ptime['sfact']
			part['redsh'] = ptime['redsh']


		return part


