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
		print('MED = ', np.median(array))
		print('AVG = ', np.average(array))

##-----
## LOAD CATALOG ROUTINES
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

