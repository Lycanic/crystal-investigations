import os
import sys
import olx
import olex_hkl
import OlexVFS
import time
import math
from io import StringIO
import ntpath

import cProfile
import pstats


#import cProfile, pstats, StringIO
#pr = cProfile.Profile()
#pr.enable()


#pr.disable()
#s = StringIO.StringIO()
#sortby = 'cumulative'
#ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
#ps.print_stats()
#print s.getvalue()

from PeriodicTable import PeriodicTable
try:
  olx.current_hklsrc
except:
  olx.current_hklsrc = None
  olx.current_hklsrc_mtime = None
  olx.current_reflections = None
  olx.current_mask = None
  olx.current_space_group = None
  olx.current_observations = None

import olex
import olex_core

import time
import cctbx_controller as cctbx_controller
from cctbx import maptbx, miller, uctbx, crystal
import iotbx
from libtbx import easy_pickle, utils

from olexFunctions import OV
from scitbx.math import distributions

from History import hist

global twin_laws_d
twin_laws_d = {}

from scitbx.math import continued_fraction
from cctbx import sgtbx, xray
from cctbx.array_family import flex
import smtbx.utils

import numpy
import itertools
import operator
import fractions
from cctbx_olex_adapter import OlexCctbxAdapter


class twin_rules:
  def __init__(self,space,twin_axis,hkl_rotation,angle,fom,threshold,rbasf=[0,0.5,0]): #changed default to bad rbasf
    self.space=space
    self.twin_axis=twin_axis
    self.hkl_rotation=hkl_rotation
    self.angle=angle
    self.fom=fom
    self.rbasf=rbasf
    self.threshold=threshold

class OlexCctbxTwinLaws(OlexCctbxAdapter):

  def __init__(self):
    #super(OlexCctbxTwinLaws, self).__init__()
    OV.registerFunction(self.run,False,'twin')
    OV.registerFunction(self.get_details,False,'twin')
    OV.registerFunction(self.run_from_gui,False,'twin')

  def run_from_gui(self, personal=0):
    print("Searching for Twin laws...")
    OV.Cursor("busy", "Searching for Twin laws...")
    personal=int(personal)
    self.run(personal) #personal
    OV.Cursor()
  def setup(self):
    global twin_laws_d
  
    OlexCctbxAdapter.__init__(self)
    if OV.GetParam('snum.refinement.use_solvent_mask'):
      txt = "Sorry, using solvent masking and twinning is not supported yet."
      print(txt)
      html = "<tr><td></td><td><b>%s</b></td></tr>" %txt
      write_twinning_result_to_disk(html)
      OV.UpdateHtml()
      return
    if "_twin" in OV.HKLSrc():
      print("It looks like your hklf file is already an hklf 5 format file") #this just checks for _twin and not for it actually being hklf5
      return
    self.fofc=self.get_fo_sq_fc()
  
    print(">> Searching for Twin Laws <<")
    use_image = ""
    self.twin_laws_d = {}
    law_txt = ""
    l = 0
    self.twin_law_gui_txt = ""
    r_list=[]
    self.filename = olx.FileName()
  
    #class-wide variable initialisation
    self.hklfile=OV.HKLSrc().rsplit('\\',1)[-1].rstrip(".hkl")
    self.model=self.xray_structure()
  
    hkl_sel_num=50
    hkl,f_calc,f_obs,f_uncertainty,leastSquare,weights=self.get_differences()
  
    self.all_hkl=hkl
    fo2 = self.reflections.f_sq_obs_filtered
    symmetry_ops=fo2.crystal_symmetry().space_group().all_ops()    
    self.symmetry_matrices=[]
    for matrix in symmetry_ops:
      self.symmetry_matrices+=[numpy.reshape(matrix.as_double_array()[0:9], (3,3))]    
    self.hkl_symmetries=hkl
    for symm in self.symmetry_matrices:
      self.hkl_symmetries=numpy.append(self.hkl_symmetries,numpy.dot(symm,hkl.T).T,axis=0)
    self.all_f_calc=f_calc
    self.set_grid()
    self.all_f_obs=f_obs
    self.weights=weights
    self.all_sigma=f_uncertainty
    self.all_leastSquare=leastSquare
    self.metrical=self.getMetrical()
    self.metrical_inv=numpy.linalg.inv(self.metrical)
    self.orthogonalization=numpy.reshape(numpy.array(self.model.unit_cell().orthogonalization_matrix()),(3,3))
    self.orthogonalization_inv=numpy.linalg.inv(self.orthogonalization)
    self.recip_orth=self.orthogonalization_inv.T
    self.recip_orth_inv=numpy.linalg.inv(self.recip_orth)
    self.overlap_threshold=OV.GetParam('snum.twinning.olex2.overlap_threshold')
    self.number_laws=OV.GetParam('snum.twinning.olex2.number_laws')
    self.cartesian_threshold=self.overlap_threshold #the threshold is provided in *cartesian* reciprocal space, to match with rotax's suggested 0.002AA
    self.cartesian_threshold_multiplier=numpy.min(numpy.linalg.norm(self.recip_orth,axis=0))
    self.relative_threshold=self.cartesian_threshold/self.cartesian_threshold_multiplier
    self.rounded_is_closest=self.get_threshold_for_rounding()
    
    rank=numpy.argsort(leastSquare)[::-1]
    hkl_sel = numpy.copy(hkl[rank[:hkl_sel_num],:])
    self.bad_fc=numpy.copy(f_calc[rank[:hkl_sel_num]])
    self.bad_fo=numpy.copy(f_obs[rank[:hkl_sel_num]])
    self.bad_hkl=hkl_sel
    self.bad_weights=numpy.copy(weights[rank[:hkl_sel_num]])
  def update_threshold(self,threshold):
    self.cartesian_threshold=threshold
    self.relative_threshold=self.cartesian_threshold/self.cartesian_threshold_multiplier
  def reset_threshold(self):
    self.cartesian_threshold=self.overlap_threshold
    self.relative_threshold=self.cartesian_threshold/self.cartesian_threshold_multiplier
  def set_grid(self):
    fcalc=self.all_f_calc
    hkl=self.all_hkl
    symm=self.hkl_symmetries
    maximums=numpy.max(abs(symm),axis=0)
    self.hkl_grid=numpy.zeros((2*int(numpy.ceil(maximums[0])),2*int(numpy.ceil(maximums[1])),2*int(numpy.ceil(maximums[2]))))
    for i,item in enumerate(symm): #note to self: symm is the list of full (non-asu) hkl, not the symmetry matrices!
      self.hkl_grid[int(item[0]),int(item[1]),int(item[2])]=fcalc[i%len(hkl)]
  def get_threshold_for_rounding(self):
    adjacents_to_zero=numpy.array([[0,0,1],[0,1,-1],[0,1,0],[0,1,1],[1,-1,-1],[1,-1,0],[1,-1,1],[1,0,-1],[1,0,0],[1,0,1],[1,1,-1],[1,1,0],[1,1,1]])
    transformed_adjacents=numpy.dot(self.recip_orth,adjacents_to_zero.T).T
    adjacent_distances=numpy.linalg.norm(transformed_adjacents,axis=1)
    minDist=numpy.min(adjacent_distances)
    maxDist=numpy.max(adjacent_distances)
    threshold_for_rounding=minDist/(2*maxDist)
    return threshold_for_rounding
  
  
  def get_details(self):
    self.setup()
    print("Orthog:")
    print(self.recip_orth)
    print("Bad Hkl")
    print(self.bad_hkl)
    
    
    

  def run(self, personal=0):
    
    self.setup()
    r_list=[]

    if personal==0:
      twin_laws=self.find_twin_laws()
    elif personal==1: #Angle-Axis
      twin_laws=self.find_AA_twin_laws()
    elif personal==2: #Sphere Search
      twin_laws=self.find_SS_twin_laws()
    else:
      print("Invalid parameter passed, personal=", personal)

    ordered_twins=sorted(twin_laws,key=lambda x: x.rbasf[1], reverse=False)
    ordered_twins=self.purge_duplicates(ordered_twins)
    top_twins=ordered_twins[:10]
    lawcount=0
    for i, twin_law in enumerate(top_twins):
      basf=twin_law.rbasf[0]
      r=twin_law.rbasf[1]
      r_diff=twin_law.rbasf[2]
      lawcount += 1
      name = " %s_twin%02d" % (self.hklfile, lawcount)
      filename = "%s.hkl" % (name)
      #self.make_hklf5(filename, twin_law, hkl, f_obs,f_uncertainty)

      self.make_hklf5_from_file(twin_law,new_file=filename) #ZZZZ
      self.twin_laws_d.setdefault(lawcount, {})

      r_no, basf_no, f_data, history = self.run_twin_ref_shelx(twin_law.hkl_rotation.flatten(), basf)

      try:
        float(r)
      except:
        r = 0.99
      r_list.append((r, lawcount, basf))

      self.twin_laws_d[lawcount]['BASF'] = basf
      self.twin_laws_d[lawcount]['r'] = r
      self.twin_laws_d[lawcount]['r_diff'] = r_diff
      self.twin_laws_d[lawcount]['matrix'] = twin_law.hkl_rotation.flatten()
      self.twin_laws_d[lawcount]['name'] = name
      self.twin_laws_d[lawcount]['HKLSrc'] = OV.HKLSrc()
    make_twinning_gui(self.twin_laws_d)




  def run_twin_ref_shelx(self, law, basf):
    #law_ins = ' '.join(str(i) for i in law[:9])
    #print "Testing: %s" %law_ins
    #file_path = olx.FilePath()
    #olx.Atreap("%s/notwin.ins" %file_path, b=True)
    #OV.AddIns("TWIN " + law_ins+" 2")
    #OV.AddIns("BASF %f"%basf)

    #curr_prg = OV.GetParam('snum.refinement.program')
    #curr_method = OV.GetParam('snum.refinement.method')
    #curr_cycles = OV.GetParam('snum.refinement.max_cycles')
    #OV.SetMaxCycles(5)
    #if curr_prg != 'olex2.refine':
    #  OV.set_refinement_program(curr_prg, 'CGLS')
    #OV.File("%s.ins" %self.filename)
    rFile = open(olx.FileFull(), 'r')
    f_data = rFile.readlines()
    rFile.close()
    #OV.SetParam('snum.init.skip_routine','True')

    #OV.SetParam('snum.refinement.program','olex2.refine')
    #OV.SetParam('snum.refinement.method','Gauss-Newton')

#    try:
#      from RunPrg import RunRefinementPrg
#      a = RunRefinementPrg()
#      self.R1 = a.R1
#      self.wR2 = a.wR2
#      his_file = a.his_file
#
#      OV.SetMaxCycles(curr_cycles)
#      OV.set_refinement_program(curr_prg, curr_method)
#    finally:
#      OV.SetParam('snum.init.skip_routine','False')


    #r = olx.Lst("R1")
    #olex_refinement_model = OV.GetRefinementModel(False)
    #if olex_refinement_model.has_key('twin'):
    #  basf = olex_refinement_model['twin']['BASF'][0]
    #else:
    #  basf = "n/a"

    return None, None, f_data, None

  def twinning_gui_def(self):
    if not self.twin_law_gui_txt:
      lines = ['search_for_twin_laws']
      tools = ['search_for_twin_laws_t1']
    else:
      lines = ['search_for_twin_laws', 'twin_laws']
      tools = ['search_for_twin_laws_t1', 'twin_laws']

    tbx = {"twinning":
           {"category":'tools',
            'tbx_li':lines
            }
           }

    tbx_li = {'search_for_twin_laws':{"category":'analysis',
                                      'image':'cctbx',
                                      'tools':['search_for_twin_laws_t1']
                                      },
              'twin_laws':{"category":'analysis',
                           'image':'cctbx',
                           'tools':['twin_laws']
                           }
              }

    tools = {'search_for_twin_laws_t1':{'category':'analysis',
                                        'display':"%Search for Twin Laws%",
                                        'colspan':1,
                                        'hrefs':['spy.OlexCctbxTwinLaws()']
                                        },
             'twin_laws':
             {'category':'analysis',
              'colspan':1,
              'before':self.twin_law_gui_txt,
              }
             }
    return {"tbx":tbx,"tbx_li":tbx_li,"tools":tools}

  def make_gui(self):
    """
    works out the figure of merit
    """
    #unused
    q = numpy.zeros((15))
    dcopy = numpy.zeros((numpy.shape(mdisag)[0]))
    fom = 0.0
    bc = 0

    #bad = numpy.rint(hkl)+numpy.rint(v) - hkl
    neighbours = list(itertools.combinations_with_replacement(numpy.arange(-1,2,dtype=numpy.float64),3))
    bad = numpy.zeros((len(neighbours), numpy.shape(mdisag)[0], 3))
    d1s = numpy.zeros((len(neighbours), numpy.shape(mdisag)[0]))
    mdisag_rint = numpy.rint(mdisag)
    for i, neighbour in enumerate(neighbours):
      bad[i,:,:] = numpy.dot(mdisag_rint+neighbour - mdisag, GS.T)
      d1s[i,:] = numpy.sqrt(numpy.sum(bad[i,:,:]*bad[i,:,:], axis=1))

    neighbours_min = numpy.zeros((numpy.shape(mdisag)[0], 3), dtype=numpy.float64)
    for i in range(numpy.shape(mdisag)[0]):
      neighbours_min[i,:] = bad[numpy.argmin(d1s[:,i]),i,:]

    bad_min = numpy.dot(mdisag_rint+neighbours_min - mdisag, GS.T)
    d_min = numpy.sqrt(numpy.sum(bad_min*bad_min, axis=1))

    dcopy = numpy.copy(d_min)
    dsum = numpy.sum(dcopy)

    for i, q_i in enumerate(q):
      q[i]=dsum/(float(numpy.shape(mdisag)[0])-float(bc))
      fom=q[i]
      if (fom<.002):
        break
      else:
        bc+=1
        jp = numpy.argmax(dcopy)
        dsum = dsum - dcopy[jp]
        dcopy[jp] = 0.

    fom=1000.*fom
    q=1000.0*q

    return fom


  def get_differences(self):
    """;
    Returns the hkl, fcalc, fobs, uncertainty and the signed difference
    """
    fot,fct=self.fofc
    weights = self.compute_weights(fot, fct)
    scale_factor = fot.scale_factor(fct, weights=weights)
    fot = fot.apply_scaling(factor=1/(scale_factor))
    leastSquare = (fot.data()-fct.norm().data())/fot.sigmas()
    hkl=numpy.array(fot.indices())
    fobs=numpy.array(fot.data())
    fcalc=numpy.array(fct.norm().data())
    funcertainty=numpy.array(fot.sigmas())
    leastSquare=numpy.array(leastSquare)
    weights=numpy.array(weights)

    return hkl,fcalc,fobs,funcertainty,leastSquare,weights

  def getMetrical(self):
    model = self.model
    gl=model.unit_cell().metrical_matrix()
    g = numpy.zeros((3,3))
    for i in range(3):
      g[i,i] = gl[i]
    g[0,1] = gl[3]
    g[0,2] = gl[4]
    g[1,2] = gl[5]
    g[1,0] = gl[3]
    g[2,0] = gl[4]
    g[2,1] = gl[5]
    return g

  def get_integral_twin_laws(self):
    twin_laws=[]
    with HiddenPrints():
      cctbx_twin_laws = cctbx_controller.twin_laws(self.reflections)

    for twin_law in cctbx_twin_laws:
      law = twin_law.as_double_array()[0:9]
      law=numpy.around(numpy.reshape(numpy.array(law),(3,3)),decimals=3)
      [basf, r,r_diff]=self.basf_estimate(law)
      if(basf<0.02):
        continue
      twin_laws+=[twin_rules("Integral",[],law,0,0,self.cartesian_threshold,[basf,r,r_diff])]

    return twin_laws

  def find_twin_laws(self):
    hkl_all=self.all_hkl
    f_calc=self.all_f_calc
    f_obs=self.all_f_obs
    hkl=self.bad_hkl

    number_laws=OV.GetParam('snum.twinning.olex2.number_laws', 4)

    twin_laws=[]


    model=self.model
    found_early=False

    tests=False
    #tests=True

    if tests==True:

      #This is a test function which replaces the lower ones
      twin_laws+=self.find_twin_laws_tests()
      if twin_laws:
        found_early=True
    else:
      #threshold=0.002
      #self.update_threshold(threshold)
      threshold=self.cartesian_threshold
      prev=len(twin_laws)
      twin_laws+=self.get_integral_twin_laws()
      post=len(twin_laws)
      if twin_laws:
        twin_laws=self.purge_duplicates(twin_laws)
        print ("Found " +str(post-prev)+" Integral Twins")

      if len(twin_laws)<number_laws:
        prev=len(twin_laws)
        twin_laws+=self.find_twofold_axes_sphere(hkl, threshold)
        post=len(twin_laws)
        if twin_laws:
          twin_laws=self.purge_duplicates(twin_laws)
          print (str(post-prev)+" Twin laws found via twofold sphere search")
      if len(twin_laws)<number_laws:
        prev=len(twin_laws)
        twin_laws+=self.find_twin_axes_sphere(hkl, threshold)
        post=len(twin_laws)
        if twin_laws:
          print(str(post-prev)+" Twin laws found through sphere search")
          twin_laws=self.purge_duplicates(twin_laws)
      if len(twin_laws)<number_laws:
        self.update_threshold(0.005)
        threshold=0.005
        prev=len(twin_laws)
        twin_laws+=self.find_twofold_axes_sphere(hkl, threshold)
        twin_laws+=self.find_twin_axes_sphere(hkl, threshold)
        post=len(twin_laws)
        if twin_laws:
          print(str(post-prev)+" Twin laws found through twofold & sphere search with large threshold")
          twin_laws=self.purge_duplicates(twin_laws)

      if len(twin_laws)<number_laws:
        rotation_fraction=2
        size=5
        prev=len(twin_laws)
        twin_laws+=self.find_twin_axes(threshold,size,rotation_fraction)
        post=len(twin_laws)
        if twin_laws:
          twin_laws=self.purge_duplicates(twin_laws)
          print (str(post-prev)+" Twin laws found via twofold rotation")

      if twin_laws:
        found_early=True

      if len(twin_laws)<number_laws and OV.GetParam('snum.twinning.olex2.do_long'): # 
        print ("Insufficient twin laws found via quicker search, trying basic brute-force search")
        olx.Refresh()
        rotation_fraction=24
        twin_laws+=self.find_twin_axes(threshold,size,rotation_fraction)

    if len(twin_laws)<number_laws and OV.GetParam('snum.twinning.olex2.do_long') and not found_early:
      print ("Checking Extended Brute-Force Search")
      olx.Refresh()
      rotation_fraction=24
      size=12
      twin_laws+=self.find_twin_axes(threshold,size,rotation_fraction)
      if twin_laws:
        twin_laws=self.purge_duplicates(twin_laws)
        if not found_early:
          print ("Twin found on extended search")

    twin_laws=self.purge_duplicates(twin_laws)
    if twin_laws:
      twin_laws=self.do_rounding(twin_laws)
      olex.m("html.ItemState * 0 tab* 2 tab-tools 1 logo1 1 index-tools* 1 info-title 1")
      olex.m("html.ItemState h2-tools-twinning 1")
      print ("Twin Laws Found - See the Twinning Tab")
      olx.Refresh()
    else:
      if not OV.GetParam('snum.twinning.olex2.do_long'):
        print ("No Twin Laws found. If there is evidence of twinning, try the 'extended' version or input your own parameters.")
      else:
        print ("No Twin Laws found. If there is evidence of twinning, try inputting your own parameters.")
    return twin_laws

  def find_SS_twin_laws(self):
    #Find the twin laws based only on the user's input parameters
    hkl=self.bad_hkl
    number_laws=self.number_laws

    twin_laws=[]
    twin_laws+=self.get_integral_twin_laws()


    size=OV.GetParam('snum.twinning.olex2.size_sphere')
    threshold=self.cartesian_threshold
    print("Using Spherical Search: \n Bad HKL to check: %d, \n Cutoff Threshold %.4f"%(size,threshold))
    olx.Refresh()

    twin_laws+=self.find_twofold_axes_sphere(hkl, threshold,size)
    if len(twin_laws)<number_laws:
      twin_laws+=self.find_twin_axes_sphere(hkl, threshold,size)

    if twin_laws:
      twin_laws=self.purge_duplicates(twin_laws)
      twin_laws=self.do_rounding(twin_laws)
      olex.m("html.ItemState * 0 tab* 2 tab-tools 1 logo1 1 index-tools* 1 info-title 1")
      olex.m("html.ItemState h2-tools-twinning 1")
      print ("Twin Laws Found - See the Twinning Tab")
      olx.Refresh()
    else:
      print ("No Twin Laws found. If you think there is likely twinning, try a larger number of bad hkl, or the angle-axis search, or ask for a greater number of laws to be considered.")
    return twin_laws


  def find_AA_twin_laws(self):
    #Find the twin laws based only on the user's input parameters
    hkl=self.bad_hkl

    twin_laws=[]
    twin_laws+=self.get_integral_twin_laws()

    rotation_fraction=OV.GetParam('snum.twinning.olex2.rotation_fraction')
    size=OV.GetParam('snum.twinning.olex2.max_index')
    threshold=self.cartesian_threshold
    print("Using personal values: \n Max Index: %d, \n Rotation Fractions: %d, \n Cutoff Threshold %.4f"%(size,rotation_fraction,threshold))
    olx.Refresh()

    twin_laws+=self.find_twin_axes(threshold,size,rotation_fraction)

    if twin_laws:
      twin_laws=self.purge_duplicates(twin_laws)
      twin_laws=self.do_rounding(twin_laws)
      olex.m("html.ItemState * 0 tab* 2 tab-tools 1 logo1 1 index-tools* 1 info-title 1")
      olex.m("html.ItemState h2-tools-twinning 1")
      print ("Twin Laws Found - See the Twinning Tab")
      olx.Refresh()
    else:
      print ("No Twin Laws found. If you think there is likely twinning, try increased index, rotation fraction or threshold.")

    return twin_laws



  def basf_estimate(self,twin_law):

    hkl=self.all_hkl
    Fc_sq=self.all_f_calc
    hkl_grid=self.hkl_grid

    hkl_new=numpy.dot(twin_law, hkl.T).T
    hkl_rint=self.closest_points(hkl_new)
    num_data=numpy.shape(hkl)[0]
    twin_component=numpy.zeros(num_data)


    fo2 = self.fofc[0]
    if numpy.allclose(twin_law,numpy.rint(twin_law)):
      small_indices=[range(0,hkl_new.shape[0])]
    else:  
      hkl_diff_euclidean=numpy.dot(self.recip_orth,(hkl_rint-hkl_new).T).T
      hkl_comparison=numpy.linalg.norm(hkl_diff_euclidean,axis=1)
      small_indices= numpy.where(hkl_comparison<=self.cartesian_threshold)
    
    for i in small_indices[0]:
      #rotated_hkl=hkl_new[i]
      hkl_new_i=hkl_rint[i]
      try:
        twin_component[i]=hkl_grid[int(hkl_new_i[0]),int(hkl_new_i[1]),int(hkl_new_i[2])]
      except:
        pass #I am sorry but the fact is 'if there isn't a value retrievable just say 0' is the current process, and doing the check in advance is a time waste. Time is vital here.
      #sometimes we are just missing a hkl, which gets set as 0 but sometimes we manage to go outside the range of the generated grid which is a problem

    if numpy.max(twin_component)==0:
      return [0,0.5,0] #basf 0, rmin=0.5, rdiff=0 - this mostly just indicates a bad twin law, although it should never reach here as laws without overlap should be ignored earlier.

    basf_r=self.find_basf_r(twin_component)

    return basf_r


  def basf_estimate_short(self,twin_law):

    hkl=self.bad_hkl
    Fc_sq=self.all_f_calc

    hkl_new=numpy.dot(twin_law, hkl.T).T
    num_data=numpy.shape(hkl)[0]
    twin_component=numpy.zeros(num_data)


    fo2 = self.fofc[0]
    if numpy.allclose(twin_law,numpy.rint(twin_law)):
      small_indices=[range(0,hkl_new.shape[0])]
    else:  
      hkl_diff_euclidean=numpy.dot(self.recip_orth,(self.closest_points(hkl_new)-hkl_new).T).T
      hkl_comparison=numpy.linalg.norm(hkl_diff_euclidean,axis=1)
      small_indices= numpy.where(hkl_comparison<=self.cartesian_threshold)
  
    hkl_all=self.all_hkl    
    for i in small_indices[0]:
      rotated_hkl=hkl_new[i]
      hkl_new_i=self.closest_points(rotated_hkl)
      loc=numpy.where((hkl_all[:,0]==hkl_new_i[0]) & (hkl_all[:,1]==hkl_new_i[1]) & (hkl_all[:,2]==hkl_new_i[2]))[0]
      if loc:
        twin_component[i]=Fc_sq[loc]
      else:
        loc=numpy.where((hkl_all[:,0]==-hkl_new_i[0]) & (hkl_all[:,1]==-hkl_new_i[1]) & (hkl_all[:,2]==-hkl_new_i[2]))[0]
        if loc:
          twin_component[i]=Fc_sq[loc]
        else:
          index=flex.miller_index()
          index.append(hkl_new_i.astype(int))
      
          #this is a botch to get the symmetry/anomaly correct as I don't understand them
          miller_set=miller.set(crystal_symmetry=fo2.crystal_symmetry(),
                                    indices=index,
                    anomalous_flag=fo2.anomalous_flag())
          #new_thing=miller.set(crystal_symmetry, index)
          miller_set=miller_set.map_to_asu()
      
          index=miller_set.indices() #can we get the 0th index easily?
          hkl_symmetric=numpy.array(index)[0]
          loc=numpy.where((hkl[:,0]==hkl_symmetric[0]) & (hkl[:,1]==hkl_symmetric[1]) & (hkl[:,2]==hkl_symmetric[2]))[0]
          if loc:
            twin_component[i]=Fc_sq[loc]           


    if numpy.max(twin_component)==0:
      return [0,0.5,0] #basf 0, rmin=0.5, rdiff=0 - this mostly just indicates a bad twin law, although it should never reach here as laws without overlap should be ignored earlier.

    basf_r=self.find_basf_r(twin_component,short=True)

    return basf_r




  def find_basf_r(self,twin_component,short=False):
    basf_lower=0
    basf_upper=1
    trials=0
    max_trials=100
    accuracy=0.01
    r_lower=self.r_from_basf(basf_lower,twin_component,short)
    r_upper=self.r_from_basf(basf_upper,twin_component,short)
    r_0=r_lower
    if r_upper>r_lower:
      #prioritise lower over upper basf value
      basf_minimum=basf_lower
      r_minimum=r_lower
    else:
      basf_minimum=basf_upper
      r_minimum=r_upper
    #need to establish a minimum, which could be at 0, before we can employ the golden section search
    while trials<max_trials and basf_upper-basf_lower>accuracy:
      basf_new=basf_lower+2/(3+math.sqrt(5))*(basf_upper-basf_lower)
      r_new=self.r_from_basf(basf_new,twin_component,short)
      if r_new<r_minimum:
        basf_minimum=basf_new
        r_minimum=r_new
        break
      #assumption - it will be less than one of the upper and lower bounds. If it were above both, it's a weird upper-peak which breaks our assumption of a single lower peak.
      elif r_new<=r_upper:
        #we only need to replace the 'upper' value, as the lower one will be the minimum if this is the case (it's not below the minimum and is below the upper)
        basf_upper=basf_new
        r_upper=r_new
      elif r_new<=r_lower:
        basf_lower=basf_new
        r_lower=r_new
      else:
        #print "error: basf midpoint has higher r value. Setting as highest basf"
        basf_upper=basf_new
        r_upper=r_new
      trials+=1

    #this only runs if we have a minimum not at 0 or 1
    while trials<max_trials and basf_upper-basf_lower>accuracy:
      basf_new=basf_lower+basf_upper-basf_minimum
      r_new=self.r_from_basf(basf_new,twin_component,short)
      if r_new<=r_minimum:
        if basf_minimum<basf_new:
          basf_lower=basf_minimum
          r_lower=r_minimum
          basf_minimum=basf_new
          r_minimum=r_new
        elif basf_minimum>basf_new:
          basf_upper=basf_minimum
          r_upper=r_minimum
          basf_minimum=basf_new
          r_minimum=r_new
        else:
          print("basf equivalent, breaking")
          break
      else: #r_new>r_minimum
        if basf_minimum<basf_new:
          basf_upper=basf_new
          r_upper=r_new
        elif basf_minimum>basf_new:
          basf_lower=basf_new
          r_lower=r_new
        else:
          print("basf equivalent, breaking")
          break
      trials +=1

    r_difference=r_0-r_minimum
    if r_difference<0.001:
      return [0,r_0,0]

    return [basf_minimum, r_minimum, r_difference]


  def r_from_basf(self,basf,twin_component,short=False):
    if short:
      Fc_sq=self.bad_fc
      Fo_sq=self.bad_fo
      w=self.bad_weights
    else:
      Fc_sq=self.all_f_calc
      Fo_sq=self.all_f_obs
      w=self.weights
    Fo_sq=numpy.maximum(Fo_sq,0)
    fcalc=(1-basf)*Fc_sq+basf*twin_component
    scale = numpy.dot(w*Fo_sq,fcalc)/numpy.dot(w*fcalc,fcalc)
    fcalc = scale * fcalc
    #R = numpy.sum(numpy.abs(numpy.sqrt(numpy.maximum(Fo_sq,0))-numpy.sqrt(fcalc)))/numpy.sum(numpy.sqrt(numpy.maximum(Fo_sq,0)))
    r1=numpy.sum(abs(abs(numpy.sqrt(Fo_sq))-abs(numpy.sqrt(fcalc))))/numpy.sum(abs(numpy.sqrt(Fo_sq)))
    wr2=numpy.sqrt(numpy.sum(abs(numpy.dot(w*(Fo_sq-fcalc),Fo_sq-fcalc)))/numpy.sum(abs(w*Fo_sq**2)))
    return r1
  
  
  #fcwfc=np.dot(w*fc_sq,fc_sq)
  #ktilde=np.dot(w*fo_sq,fc_sq)/fcwfc   
  #r1=np.sum(abs(abs(np.sqrt(fo_sq))-abs(np.sqrt(ktilde*fc_sq))))/np.sum(abs(np.sqrt(fo_sq)))
  #r3=np.sum(abs(np.dot(w*(fo_sq-ktilde*fc_sq),fo_sq-ktilde*fc_sq)))/np.sum(abs(w*fo_sq**2))
  #wr2=np.sqrt(r3) ##this is the wR2  


  def find_twin_axes(self, threshold,size,rotation_fraction):

    hkl=self.all_hkl
    model=self.model

    metrical_matrix=self.metrical
    metrical_inverse=self.metrical_inv
    orthogonalization_matrix=self.orthogonalization
    orthogonalization_inverse=self.orthogonalization_inv
    reciprocal_orthogonalization_matrix=self.recip_orth
    reciprocal_orthogonalization_inverse=self.recip_orth_inv
    self.update_threshold(threshold)
  

    perfect_reflection_number=math.exp(-13)
    possible_twin_laws=[]

    #angle and the sine and cosine of that, for use in the rotation formula
    base_rotation_angle=2.*math.pi/rotation_fraction

    for twin_axis in self.all_axes(size): #itertools.product(numpy.arange(-size,size+1),numpy.arange(-size,size+1),range(size+1)):
          reciprocal_law=False
          ##skip inverse axes
          #if(twin_axis[2]==0):
            #if(twin_axis[1]<0):
              #continue
            #if(twin_axis[1]==0):
              #if(twin_axis[0]<=0):
                #continue
          #if(fractions.gcd(fractions.gcd(twin_axis[0],twin_axis[1]),twin_axis[2])!=1):
            #continue

          #using the rodrigues formula to generate the matrices, the O^-1RO for that in lattice coordinates
          rotation_matrix_lattice_base=self.make_lattice_rotation(twin_axis,base_rotation_angle,orthogonalization_matrix,reciprocal_orthogonalization_matrix)
          recip_rot_lat_base=self.make_lattice_rotation(twin_axis,base_rotation_angle,reciprocal_orthogonalization_matrix,reciprocal_orthogonalization_matrix)

          new_hkl=hkl.copy()
          rec_hkl=hkl.copy()

          #this is the rotation matrix in the lattice coordinates! due to quirks and transforming a vector into cartesian space, rotating and
          #transforming back being a O^-1 R O type, we can square and more for 'bigger' angles.
          for r in numpy.arange(1,rotation_fraction): #every part of this loop relies on the previous completing the hkl rotations - this is for efficiency (about 10%)

            #reciprocal axis
            rec_hkl=numpy.dot(rec_hkl,recip_rot_lat_base.T)
            rec_hkl_displacement=self.find_fom(rec_hkl)

            #real axis
            new_hkl=numpy.dot(rotation_matrix_lattice_base,new_hkl.T).T #correct
            #if not self.sufficient_overlaps(new_hkl, threshold) and not self.sufficient_overlaps(rec_hkl, threshold) :
              #continue
            hkl_displacement=self.find_fom(new_hkl)

            if (rec_hkl_displacement<threshold and rec_hkl_displacement>perfect_reflection_number):
              reciprocal_rotation_lattice=numpy.linalg.matrix_power(recip_rot_lat_base,r)
              rbasf=self.basf_estimate(reciprocal_rotation_lattice)
              if rbasf[0]>1e-5:
                possible_twin_laws+=[twin_rules("Reciprocal",twin_axis,numpy.around(reciprocal_rotation_lattice,decimals=3),r*base_rotation_angle,rec_hkl_displacement,self.cartesian_threshold,rbasf=rbasf)]
                reciprocal_law=True

            if (hkl_displacement<threshold and hkl_displacement>perfect_reflection_number):
              rotation_matrix_lattice=numpy.linalg.matrix_power(rotation_matrix_lattice_base,r)
              rbasf=self.basf_estimate(rotation_matrix_lattice)
              if rbasf[0]>1e-5:
                possible_twin_laws+=[twin_rules("Direct",twin_axis,numpy.around(rotation_matrix_lattice,decimals=3),r*base_rotation_angle,hkl_displacement,self.cartesian_threshold,rbasf=rbasf)]
    self.reset_threshold()
    return possible_twin_laws

  def make_lattice_rotation(self, axis, angle,axis_orthog, reciprocal_orthogonalization_matrix):
    cosine_value=math.cos(angle)
    sine_value=math.sin(angle)
    reciprocal_orthogonalization_inverse=numpy.linalg.inv(reciprocal_orthogonalization_matrix)

    axis_cartesian=numpy.dot(axis_orthog,axis)
    axis_unit_cartesian=axis_cartesian/self.size_of_3d_vector(axis_cartesian)
    cross_product_matrix=numpy.array([[0,-axis_unit_cartesian[2],axis_unit_cartesian[1]],[axis_unit_cartesian[2],0,-axis_unit_cartesian[0]],[-axis_unit_cartesian[1],axis_unit_cartesian[0],0]],dtype=float)
    matrix=numpy.eye(3)+sine_value*cross_product_matrix+(1-cosine_value)*numpy.linalg.matrix_power(cross_product_matrix,2)
    matrix_lattice=numpy.dot(reciprocal_orthogonalization_inverse,numpy.dot(matrix,reciprocal_orthogonalization_matrix))

    return matrix_lattice


  def make_real_lattice_rotation(self, axis, angle): #believe this is actually the euclidean rotation - thus why I called it 'real'
    cosine_value=math.cos(angle)
    sine_value=math.sin(angle)

    cross_product_matrix=numpy.array([[0,-axis[2],axis[1]],[axis[2],0,-axis[0]],[-axis[1],axis[0],0]],dtype=float)
    matrix=numpy.eye(3)+sine_value*cross_product_matrix+(1-cosine_value)*numpy.linalg.matrix_power(cross_product_matrix,2)

    return matrix

  def find_fom(self, falsehkl):
    metrical=self.metrical_inv
    if all(numpy.sum(metrical, axis=1)>=0): #this is a loose condition which dictates rounding to the nearest integer will also give the closest reciprocal lattice point.
      displacement = (falsehkl+0.5)%1-0.5
      distances=numpy.sqrt(numpy.multiply(displacement,numpy.dot(metrical,displacement.T).T).sum(1)) #      distances=numpy.linalg.norm(numpy.dot(orthogonalization.T,displacement.T),axis=0) is alternative but slower
    else:
      hkl_dropped=numpy.floor(falsehkl)
      hkl_remainder=falsehkl%1-1 #where it would round to the closest, 1,1,1 will provide that distance!
      adjacents=numpy.array([[0,0,0],[1,0,0],[0,1,0],[1,1,0],[0,0,1],[1,0,1],[0,1,1],[1,1,1]])
      distances=numpy.full(numpy.shape(falsehkl)[0],1000)
      for i in numpy.arange(0,8):
        displacement=hkl_remainder+adjacents[i]
        #hkl_adjacent=hkl_dropped+adjacents[i]
        #displacement=falsehkl-hkl_adjacent
        size=numpy.sqrt(numpy.multiply(displacement,numpy.dot(metrical,displacement.T).T).sum(1))
        distances=numpy.minimum(distances,size)
    sortedDist=numpy.sort(distances)
    length=len(distances)
    filtered=sortedDist[:-int(length/2)]
    return numpy.average(filtered)
  def closest_points(self,hkl): # returns the closest points to hkl, bearing in mind that if the overlap threshold is sufficiently small then all points within that threshold will be mapped to their rounded point
    if self.overlap_threshold<=self.rounded_is_closest:
      return numpy.rint(hkl)
    else:
      hkl_dropped=numpy.floor(hkl) # hkl = dropped+remainder. Possibilities: dropped + delta. Distance = delta-remainder?
      hkl_remainder=hkl%1 #provides a decimal in 0 to 1
      adjacents=numpy.array([[0,0,0],[1,0,0],[0,1,0],[1,1,0],[0,0,1],[1,0,1],[0,1,1],[1,1,1]])
      distances=numpy.full(numpy.shape(hkl)[0],1000)
      closest=numpy.full(numpy.shape(hkl),0)
      for i in numpy.arange(0,8):
        displacement=-hkl_remainder+adjacents[i]
        size=numpy.sqrt(numpy.multiply(displacement,numpy.dot(self.metrical_inv,displacement.T).T).sum(1))
        distances=numpy.minimum(distances,size)      
        bestLocations=numpy.where(distances==size)
        closest[bestLocations]=hkl_dropped[bestLocations]+adjacents[i]
      return closest
  def find_number_overlaps(self, falsehkl, orthogonalization, threshold):
    metrical=self.metrical_inv
    displacement = (falsehkl+0.5)%1-0.5
    threshold2=threshold**2
    overlaps=0
    size=numpy.sqrt(numpy.multiply(displacement,numpy.dot(metrical,displacement.T).T).sum(1))#size=numpy.dot(displacement,numpy.dot(metrical,displacement.T))
    for line in displacement:
      size=numpy.dot(line,numpy.dot(metrical,line))
      if size<threshold2:
        overlaps+=1
    return overlaps


  def sufficient_overlaps(self,hkl,threshold):
    #returns 0 if the number of overlaps drops below a certain amount (nominally 35 of 50), else returns the number of overlaps
    #the threshold is in relative coordinates currently
    metrical=self.metrical_inv
    displacement = (hkl+0.5)%1-0.5
    threshold2=threshold**2
    overlaps=50
    for line in displacement:
      size=numpy.dot(line,numpy.dot(metrical,line))
      if size>threshold2:
        overlaps-=1
      if overlaps<35:
        return 0
    return overlaps

  def make_hklf5(self, filename, twin_law_full, hkl, fo,sigmas):
    """
    convert hklf4 file to hklf5 file
    """
    #Pulled directly from Pascal's, then edited to generate its own second component hkl, as some files don't have avlues for all


    twin_law=twin_law_full.hkl_rotation
    basf=twin_law_full.rbasf[0]
    cell=self.cell
    hklf = open(filename,'w')
    hkl_new = numpy.dot(twin_law, hkl.T).T
    rounded_hkl_new=self.closest_points(hkl_new).astype(int)

    scale = 99999/numpy.max(fo) # keep Fo in the scale of F8
    threshold=self.cartesian_threshold
    check1=True

    for i, hkl_new_i in enumerate(hkl_new):
      if(numpy.any(numpy.abs(self.closest_points(hkl_new_i)-hkl_new_i)>0.1)):
        hklf.write("%4d%4d%4d%8d%8d%4d\n"%(hkl[i,0],hkl[i,1], hkl[i,2], fo[i]*scale, sigmas[i]*scale, 1))
        continue
      # looking where is the overlap
      #diff = numpy.sum(numpy.abs(hkl_new[i] - rounded_hkl_new[i]))
      diff = numpy.linalg.norm(numpy.dot(self.recip_orth(hkl_new[i] - rounded_hkl_new[i])))

      if(diff<threshold): #threshold prior default is 0.1
        # adding contribution from overlap
        hklf.write("%4d%4d%4d%8d%8d%4d\n"%(rounded_hkl_new[i,0],rounded_hkl_new[i,1], rounded_hkl_new[i,2], fo[i]*scale, sigmas[i]*scale, -2))
        hklf.write("%4d%4d%4d%8d%8d%4d\n"%(hkl[i,0],hkl[i,1], hkl[i,2], fo[i]*scale, sigmas[i]*scale, 1))
      else:
        hklf.write("%4d%4d%4d%8d%8d%4d\n"%(hkl[i,0],hkl[i,1], hkl[i,2], fo[i]*scale, sigmas[i]*scale, 1))

    hklf.write("%4d%4d%4d\n"%(0,0,0))
    hklf.write("CELL %s %s %s %s %s %s\n"%(cell[0],cell[1],cell[2],cell[3],cell[4],cell[5]))
    hklf.write("REM TWIN %s\n"%format_twin_string_from_law(twin_law))
    hklf.write("BASF %6.4f\n"%basf)
    hklf.write("HKLF 5")
    hklf.close()

  def make_hklf5_from_file(self, twin_law_full, filename=None,new_file=None):
    """
    convert hklf4 file to hklf5 file using the file as a base
    """
    if not filename:
      filename=self.hklfile+'.hkl' #OV.HKLSrc().rsplit('\\',1)[-1]
    if not new_file:
      new_file=self.hklfile+'twinX.hkl'
    metrical=self.getMetrical()
    metrical_inv=numpy.linalg.inv(metrical)
    twin_law=twin_law_full.hkl_rotation
    hkl_file_old=list(olex_hkl.Read(filename))
    hkl_file_new=[]
    for line in hkl_file_old:
      line_as_list=list(line)
      [h,k,l,fo,sig,batch,nul2]=line_as_list
      hkl=[h,k,l]
      hkl_new = numpy.dot(twin_law, hkl)
      hkl_new_rounded=self.closest_points(hkl_new).astype(int)
      #difference=hkl_new_rounded-hkl_new
      #distance=numpy.sqrt(numpy.dot(difference,numpy.dot(metrical_inv,difference)))
      #if distance<threshold:
      #if not (numpy.any(numpy.abs(self.closest_points(hkl_new)-hkl_new)>0.1)):
        #diff = numpy.sum(numpy.abs(self.closest_points(hkl_new)-hkl_new))
      diff = numpy.linalg.norm(numpy.dot(self.recip_orth,(hkl_new_rounded-hkl_new)))
      if diff<twin_law_full.threshold:
        hkl_file_new+=[(hkl_new_rounded[0],hkl_new_rounded[1],hkl_new_rounded[2],fo,sig,-2)]
      hkl_file_new+=[(h,k,l,fo,sig,1)]
    try:
      olex_hkl.Write(new_file,hkl_file_new)
    except Exception:
      pass

    basf=twin_law_full.rbasf[0]
    cell=self.cell
    hklf = open(new_file,'a')

    hklf.write("%4d%4d%4d\n"%(0,0,0))
    hklf.write("CELL %s %s %s %s %s %s\n"%(cell[0],cell[1],cell[2],cell[3],cell[4],cell[5]))
    hklf.write("REM TWIN %s\n"%format_twin_string_from_law(twin_law))
    hklf.write("BASF %6.4f\n"%basf)
    hklf.write("HKLF 5")
    hklf.close()

  def extend_hklf5(self, twin_law_full,new_component, filename=None,new_file=None, threshold=0.002, component_base=1):
    """
    adds an extra twin law to a hklf 5 or 4 file. The twin_law_full rotates the component_base component, and will be checked against that component only. The new component has integer new_component, which should differ from any other batch numbers
    """
    if not filename:
      filename=self.hklfile+'.hkl'
    if not new_file:
      new_file=self.hklfile+'twinX.hkl'
    metrical=self.getMetrical()
    metrical_inv=numpy.linalg.inv(metrical)
    twin_law=twin_law_full.hkl_rotation
    hkl_file_old=list(olex_hkl.Read(filename))
    hkl_file_new=[]
    for line in hkl_file_old:
      line_as_list=list(line)
      [h,k,l,fo,sig,batch,nul2]=line_as_list
      hkl=[h,k,l]
      hkl_new = numpy.dot(twin_law, hkl)
      hkl_new_rounded=self.closest_points(hkl_new).astype(int)
      #difference=hkl_new_rounded-hkl_new
      #distance=numpy.sqrt(numpy.dot(difference,numpy.dot(metrical_inv,difference)))
      #if distance<threshold:
      if abs(batch)==component_base:
        if not (numpy.any(numpy.abs(self.closest_points(hkl_new)-hkl_new)>0.1)):
          #diff = numpy.sum(numpy.abs(self.closest_points(hkl_new)-hkl_new))
          diff = numpy.linalg.norm(numpy.dot(self.recip_orth(self.closest_points(hkl_new)-hkl_new)))
          if diff<self.cartesian_threshold:
            hkl_file_new+=[(hkl_new_rounded[0],hkl_new_rounded[1],hkl_new_rounded[2],fo,sig,-new_component)]
        hkl_file_new+=[(h,k,l,fo,sig,1)]
    try:
      olex_hkl.Write(new_file,hkl_file_new)
    except Exception:
      pass

    basf=twin_law_full.rbasf[0]
    cell=self.cell
    hklf = open(new_file,'a')

    hklf.write("%4d%4d%4d\n"%(0,0,0))
    hklf.write("CELL %s %s %s %s %s %s\n"%(cell[0],cell[1],cell[2],cell[3],cell[4],cell[5]))
    hklf.write("REM TWIN %s\n"%format_twin_string_from_law(twin_law))
    hklf.write("BASF %6.4f\n"%basf)
    hklf.write("HKLF 5")
    hklf.close()






  def do_rounding(self,twin_laws):
    for twin_law in twin_laws:
      matrix=twin_law.hkl_rotation
      for i in numpy.arange(0,3):
        for j in numpy.arange(0,3):
          item=matrix[i,j]
          if (abs(1.0-item)<0.01):
            matrix[i,j]=1.0
          elif (abs(-1.0-item)<0.01):
            matrix[i,j]=-1.0
          elif (abs(0.0-item)<0.01):
            matrix[i,j]=0.0
          elif (abs(0.5-item)<0.01):
            matrix[i,j]=0.5
      twin_law.hkl_rotation=matrix
    return twin_laws

  def purge_duplicates(self, twin_laws):
    twin_laws=self.check_basf(twin_laws)
    twin_laws=sorted(twin_laws,key=lambda x: x.rbasf[1], reverse=False)

    non_duplicate_twin_laws=[]
    fo2 = self.reflections.f_sq_obs_filtered
    #symmetry_ops=fo2.crystal_symmetry().space_group().all_ops()#list of symmetries to be retrieved with .as_double_array() and mapped to actual 3x3 matrices.
    #However, they have 3 more numbers than expected - the last 3 are the +- for translations I think. Sadly, due to inversion symmetry, there are twice as many matrices as actual symmetry equivalents we need to care about.
    #symmetry_matrices=[]
    #for matrix in symmetry_ops:
      #symmetry_matrices+=[numpy.reshape(matrix.as_double_array()[0:9], (3,3))]
    symmetry_matrices=self.symmetry_matrices
    #assumption - only the first half are unique, the rest represent inversion symmetries, which are inherantly discarded by only taking positive axes.
    #additionally, the first is the identity, so is also discarded
    #num_matrix=len(symmetry_matrices)
    #symmetry_matrices=symmetry_matrices[0:int(num_matrix/2)] #changed from 1 to 0 so it purges pure duplicates

    #I have now removed the removal of inversion, seeing as I'm using a different method it does not apply
    for i, twin_law in enumerate(twin_laws):
      if i==0:
        non_duplicate_twin_laws+=[twin_law]
        continue
      duplicate=False
      for matrix in symmetry_matrices:
        rotated_twin=numpy.dot(matrix,twin_law.hkl_rotation)
        for j in range(0,i):
          if numpy.allclose(rotated_twin, twin_laws[j].hkl_rotation,atol=5e-3): #same up to 2dp.
            duplicate=True
          if duplicate:
            break;
        if duplicate:
          break;
      if not duplicate:
        non_duplicate_twin_laws+=[twin_law]


    return non_duplicate_twin_laws

  def all_axes(self, max_index):
    fo2 = self.reflections.f_sq_obs_filtered
    indexes=flex.miller_index()
    for twin_axis in itertools.product(numpy.arange(-max_index,max_index+1),numpy.arange(-max_index,max_index+1),list(range(max_index+1))):
      indexes.append(twin_axis)

    miller_set=miller.set(crystal_symmetry=fo2.crystal_symmetry(),
          indices=indexes,
          anomalous_flag=fo2.anomalous_flag())
    miller_set=miller_set.map_to_asu()
    asu_indexes=numpy.array(miller_set.indices())
    asu_unique=numpy.unique(asu_indexes,axis=0)
    for i in range(asu_unique.shape[0]):
      if all(asu_unique[i]==[0,0,0]) or fractions.gcd(asu_unique[i][0],fractions.gcd(asu_unique[i][1],asu_unique[i][2]))!=1:
        continue
      yield asu_unique[i]

  def find_twin_laws_tests(self):
    filename=OV.FileName()+'_test_records.txt'
    open(filename,'w').close()
    hkl=self.bad_hkl
    hkl_all=self.all_hkl
    hkl_symmetries=self.hkl_symmetries
    hkl_grid=self.hkl_grid

    twin_laws=[]
    import cProfile, pstats, io
    sortby = 'cumulative'
    
    size_limit=10
    RF=2
    threshold_finding=0.002
    self.update_threshold(threshold_finding)



    pr = cProfile.Profile()
    pr.enable()

    laws=self.get_integral_twin_laws()

    pr.disable()
    s = StringIO()
    ps = pstats.Stats(pr, stream=s).sort_stats(sortby)

    ps.print_stats(100)

    testFile=open(filename,'a')
    testFile.write("CCtbx Integral Twins\n")
    testFile.write("Laws found: "+str(len(laws))+"\n")
    self.tests_writer(laws,testFile)
    testFile.write(s.getvalue())
    testFile.close()

    twin_laws+=laws


    pr = cProfile.Profile()
    pr.enable()

    laws=self.find_twofold_axes_sphere(hkl, 0.002)#0.01)

    pr.disable()
    s = io.StringIO()
    ps = pstats.Stats(pr, stream=s).sort_stats(sortby)

    ps.print_stats(100)

    testFile=open(filename,'a')
    testFile.write("Twofold Sphere Search\n")
    testFile.write("threshold: "+str(0.002)+"\n")
    testFile.write("Laws found: "+str(len(laws))+"\n")
    self.tests_writer(laws,testFile)
    testFile.write(s.getvalue())
    testFile.close()

    twin_laws+=laws




    pr = cProfile.Profile()
    pr.enable()

    laws=self.find_twin_axes(threshold_finding, size_limit, RF)

    pr.disable()
    s = io.StringIO()
    ps = pstats.Stats(pr, stream=s).sort_stats(sortby)

    ps.print_stats(100)

    testFile=open(filename,'a')
    testFile.write("Angle-Axis Search\n")
    testFile.write("Threshold: "+str(threshold_finding)+", size: "+str(size_limit)+", Rotation Fraction: "+str(RF)+"\n")
    testFile.write("Laws found: "+str(len(laws))+"\n")
    self.tests_writer(laws,testFile)
    testFile.write(s.getvalue())
    testFile.close()
    #print s.getvalue()

    twin_laws+=laws


    ##threshold_finding=0.1
    ##self.update_threshold(threshold_finding)
    #size_limit=1


    #pr = cProfile.Profile()
    #pr.enable()

    #laws=self.find_twin_axes(threshold_finding, size_limit, RF)

    #pr.disable()
    #s = io.StringIO()
    #ps = pstats.Stats(pr, stream=s).sort_stats(sortby)

    #ps.print_stats(100)

    #testFile=open(filename,'a')
    #testFile.write("Angle-Axis Search\n")
    #testFile.write("Threshold: "+str(threshold_finding)+", size: "+str(size_limit)+", Rotation Fraction: "+str(RF)+"\n")
    #testFile.write("Laws found: "+str(len(laws))+"\n")
    #self.tests_writer(laws,testFile)
    #testFile.write(s.getvalue())
    #testFile.close()
    ##print s.getvalue()

    #twin_laws+=laws


    threshold_finding=0.002
    self.update_threshold(threshold_finding)

    pr = cProfile.Profile()
    pr.enable()

    laws=self.find_twin_axes_sphere(hkl, threshold_finding)

    pr.disable()
    s = io.StringIO()
    ps = pstats.Stats(pr, stream=s).sort_stats(sortby)

    ps.print_stats(100)

    testFile=open(filename,'a')
    testFile.write("Sphere Search\n")
    testFile.write("threshold: "+str(threshold_finding)+"\n")
    testFile.write("Laws found: "+str(len(laws))+"\n")
    self.tests_writer(laws,testFile)
    testFile.write(s.getvalue())
    testFile.close()

    twin_laws+=laws



    #threshold_finding=0.1
    #self.update_threshold(threshold_finding)

    #pr = cProfile.Profile()
    #pr.enable()

    #laws=self.find_twin_axes_sphere(hkl, threshold_finding)

    #pr.disable()
    #s = io.StringIO()
    #ps = pstats.Stats(pr, stream=s).sort_stats(sortby)

    #ps.print_stats(100)

    #testFile=open(filename,'a')
    #testFile.write("Sphere Search\n")
    #testFile.write("threshold: "+str(threshold_finding)+"\n")
    #testFile.write("Laws found: "+str(len(laws))+"\n")
    #self.tests_writer(laws,testFile)
    #testFile.write(s.getvalue())
    #testFile.close()

    #twin_laws+=laws


    hkl=self.bad_hkl
    rotated_hkl=[]
      
    for law in twin_laws:
      rotated_hkl+=[numpy.dot(hkl_all,law.hkl_rotation.T)]

    pr = cProfile.Profile()
    pr.enable()
    
    for i in range(10000):
    
      rotated_hkl=[]
      
      for law in twin_laws:
        rotated_hkl+=[numpy.dot(hkl,law.hkl_rotation.T)]
      

    pr.disable()
    s = io.StringIO()
    ps = pstats.Stats(pr, stream=s).sort_stats(sortby)

    ps.print_stats(100)  
    testFile=open(filename,'a')
    testFile.write("rotation of bad hkls\n")
    testFile.write(s.getvalue())  
    testFile.close()

    pr = cProfile.Profile()
    pr.enable()
    
    for i in range(10000):
    
      rotated_hkl=[]
      
      for law in twin_laws:
        rotated_hkl+=[numpy.dot(hkl_all,law.hkl_rotation.T)]
      

    pr.disable()
    s = io.StringIO()
    ps = pstats.Stats(pr, stream=s).sort_stats(sortby)

    ps.print_stats(100)  
    testFile=open(filename,'a')
    testFile.write("rotation of all hkls\n")
    testFile.write(s.getvalue())  
    testFile.close()


    pr = cProfile.Profile()
    pr.enable()

    for i in range(100):
      for law in twin_laws:
        self.basf_estimate(law.hkl_rotation)
        

    pr.disable()
    s = io.StringIO()
    ps = pstats.Stats(pr, stream=s).sort_stats(sortby)

    ps.print_stats(100)  
    testFile=open(filename,'a')
    testFile.write("basf tests only\n")
    testFile.write(s.getvalue())  
    testFile.close()



    #pr = cProfile.Profile()
    #pr.enable()

    #for i in range(100):
      #for law in twin_laws:
        #self.basf_estimate_short(law.hkl_rotation)
      

    #pr.disable()
    #s = io.StringIO()
    #ps = pstats.Stats(pr, stream=s).sort_stats(sortby)

    #ps.print_stats(100)  
    #testFile=open(filename,'a')
    #testFile.write("short basf tests only, state 1 and 2\n")
    #testFile.write(s.getvalue())  
    #testFile.close()



    pr = cProfile.Profile()
    pr.enable()

    for i in range(10000):
      for hklx in rotated_hkl:
        self.find_fom(hklx)
      

    pr.disable()
    s = io.StringIO()
    ps = pstats.Stats(pr, stream=s).sort_stats(sortby)

    ps.print_stats(100)  
    testFile=open(filename,'a')
    testFile.write("fom tests only\n")
    testFile.write(s.getvalue()) 
    testFile.close() 


    #pr = cProfile.Profile()
    #pr.enable()

    #for i in range(10):
      #for hklx in rotated_hkl:
        #for hkli in hklx:
          #place=numpy.min(numpy.append(numpy.where((hkl_all[:,0]==hkli[0]) & (hkl_all[:,1]==hkli[1]) & (hkl_all[:,2]==hkli[2]))[0],numpy.inf))

    #pr.disable()
    #s = io.StringIO()
    #ps = pstats.Stats(pr, stream=s).sort_stats(sortby)

    #ps.print_stats(100)  
    #testFile=open(filename,'a')
    #testFile.write("One hkl-list where search- length = "+str(len(hkl_all))+"\n")
    #testFile.write(s.getvalue()) 
    #testFile.close() 


    #pr = cProfile.Profile()
    #pr.enable()

    #for i in range(10):
      #for hklx in rotated_hkl:
        #for hkli in hklx:
          #place=numpy.min(numpy.append(numpy.where((hkl_symmetries[:,0]==hkli[0]) & (hkl_symmetries[:,1]==hkli[1]) & (hkl_symmetries[:,2]==hkli[2]))[0],numpy.inf))
      

    #pr.disable()
    #s = io.StringIO()
    #ps = pstats.Stats(pr, stream=s).sort_stats(sortby)

    #ps.print_stats(100)  
    #testFile=open(filename,'a')
    #testFile.write("symmetry-extended where - length = "+str(len(hkl_symmetries))+"\n")
    #testFile.write(s.getvalue()) 
    #testFile.close() 


    #def retrieve_point(hkli):
      #fcalc=hkl_grid[int(hkli[0]),int(hkli[1]),int(hkli[2])]
      #return fcalc
    #pr = cProfile.Profile()
    #pr.enable()

    #for i in range(1000):
      #for hklx in rotated_hkl:
        #for hkli in hklx:
          #fcalc=retrieve_point(hkli)
          

    #pr.disable()
    #s = io.StringIO()
    #ps = pstats.Stats(pr, stream=s).sort_stats(sortby)

    #ps.print_stats(100)  
    #testFile=open(filename,'a')
    #testFile.write("grid-based fcalc getting - length = "+str(hkl_grid.shape)+"\n")
    #testFile.write(s.getvalue()) 
    #testFile.close() 


    pr = cProfile.Profile()
    pr.enable()

    for i in range(10):
      self.setup()
          

    pr.disable()
    s = io.StringIO()
    ps = pstats.Stats(pr, stream=s).sort_stats(sortby)

    ps.print_stats(100)  
    testFile=open(filename,'a')
    testFile.write("setup time")
    testFile.write(s.getvalue()) 
    testFile.close() 

    #pr = cProfile.Profile()
    #pr.enable()

    #for i in range(1000):
      #self.set_grid()

    #pr.disable()
    #s = io.StringIO()
    #ps = pstats.Stats(pr, stream=s).sort_stats(sortby)

    #ps.print_stats(100)  
    #testFile=open(filename,'a')
    #testFile.write("grid creation test")
    #testFile.write(s.getvalue()) 
    #testFile.close() 

    #pr = cProfile.Profile()
    #pr.enable()

    #for i in range(10000):
      #for hkl in rotated_hkl:
        #self.sufficient_overlaps(hkl,0.01)
      

    #pr.disable()
    #s = io.StringIO()
    #ps = pstats.Stats(pr, stream=s).sort_stats(sortby)

    #ps.print_stats(100)  
    #testFile=open(filename,'a')
    #testFile.write("enough_overlaps tests only\n")
    #testFile.write(s.getvalue())  
    #testFile.close()

    #pr = cProfile.Profile()
    #pr.enable()

    #for i in range(1000):
      #for hkl in rotated_hkl:
        #self.size_of_3d_vector(hkl)
        #numpy.linalg.norm(hkl)
      
      

    #pr.disable()
    #s = io.StringIO()
    #ps = pstats.Stats(pr, stream=s).sort_stats(sortby)

    #ps.print_stats(100)  
    #testFile=open(filename,'a')
    #testFile.write("Numpy versus size\n")
    #testFile.write(s.getvalue())  
    #testFile.close()




    #RF=24

    #pr = cProfile.Profile()
    #pr.enable()

    #laws=self.find_twin_axes(threshold_finding, size_limit, RF)

    #pr.disable()
    #s = io.StringIO()
    #ps = pstats.Stats(pr, stream=s).sort_stats(sortby)

    #ps.print_stats(100)

    #testFile=open(filename,'a')
    #testFile.write("Angle-Axis Search\n")
    #testFile.write("Threshold: "+str(threshold_finding)+", size: "+str(size_limit)+", Rotation Fraction: "+str(RF)+"\n")
    #testFile.write("Laws found: "+str(len(laws))+"\n")
    #self.tests_writer(laws,testFile)
    #testFile.write(s.getvalue())
    #testFile.close()

    #twin_laws+=laws


    #threshold_finding=0.1
    #self.update_threshold(threshold_finding)

    #pr = cProfile.Profile()
    #pr.enable()

    #laws=self.find_twin_axes_sphere(hkl, threshold_finding)

    #pr.disable()
    #s = io.StringIO()
    #ps = pstats.Stats(pr, stream=s).sort_stats(sortby)

    #ps.print_stats(100)

    #testFile=open(filename,'a')
    #testFile.write("Sphere Search\n")
    #testFile.write("threshold: "+str(threshold_finding)+"\n")
    #testFile.write("Laws found: "+str(len(laws))+"\n")
    #self.tests_writer(laws,testFile)
    #testFile.write(s.getvalue())
    #testFile.close()

    #twin_laws+=laws




    #RF=24
    #size_limit=12
    ##threshold_finding=0.01
    ##self.update_threshold(threshold_finding)

    #pr = cProfile.Profile()
    #pr.enable()

    #laws=self.find_twin_axes(threshold_finding, size_limit, RF)

    #pr.disable()
    #s = io.StringIO()
    #ps = pstats.Stats(pr, stream=s).sort_stats(sortby)

    #ps.print_stats(100)

    #testFile=open(filename,'a')
    #testFile.write("Angle-Axis Search\n")
    #testFile.write("Threshold: "+str(threshold_finding)+", size: "+str(size_limit)+", Rotation Fraction: "+str(RF)+"\n")
    #testFile.write("Laws found: "+str(len(laws))+"\n")
    #self.tests_writer(laws,testFile)
    #testFile.write(s.getvalue())
    #testFile.close()





    #twin_laws+=laws
    

    threshold_finding=0.005
    self.update_threshold(threshold_finding)



    pr = cProfile.Profile()
    pr.enable()

    laws=self.get_integral_twin_laws()

    pr.disable()
    s = StringIO()
    ps = pstats.Stats(pr, stream=s).sort_stats(sortby)

    ps.print_stats(100)

    testFile=open(filename,'a')
    testFile.write("CCtbx Integral Twins\n")
    testFile.write("Laws found: "+str(len(laws))+"\n")
    self.tests_writer(laws,testFile)
    testFile.write(s.getvalue())
    testFile.close()

    twin_laws+=laws


    pr = cProfile.Profile()
    pr.enable()

    laws=self.find_twofold_axes_sphere(hkl, threshold_finding)

    pr.disable()
    s = io.StringIO()
    ps = pstats.Stats(pr, stream=s).sort_stats(sortby)

    ps.print_stats(100)

    testFile=open(filename,'a')
    testFile.write("Twofold Sphere Search\n")
    testFile.write("threshold: "+str(threshold_finding)+"\n")
    testFile.write("Laws found: "+str(len(laws))+"\n")
    self.tests_writer(laws,testFile)
    testFile.write(s.getvalue())
    testFile.close()

    twin_laws+=laws


    RF=24

    pr = cProfile.Profile()
    pr.enable()

    laws=self.find_twin_axes(threshold_finding, size_limit, RF)

    pr.disable()
    s = io.StringIO()
    ps = pstats.Stats(pr, stream=s).sort_stats(sortby)

    ps.print_stats(100)

    testFile=open(filename,'a')
    testFile.write("Angle-Axis Search\n")
    testFile.write("Threshold: "+str(threshold_finding)+", size: "+str(size_limit)+", Rotation Fraction: "+str(RF)+"\n")
    testFile.write("Laws found: "+str(len(laws))+"\n")
    self.tests_writer(laws,testFile)
    testFile.write(s.getvalue())
    testFile.close()
    #print s.getvalue()

    twin_laws+=laws



    pr = cProfile.Profile()
    pr.enable()

    laws=self.find_twin_axes_sphere(hkl, threshold_finding)

    pr.disable()
    s = io.StringIO()
    ps = pstats.Stats(pr, stream=s).sort_stats(sortby)

    ps.print_stats(100)

    testFile=open(filename,'a')
    testFile.write("Sphere Search\n")
    testFile.write("threshold: "+str(threshold_finding)+"\n")
    testFile.write("Laws found: "+str(len(laws))+"\n")
    self.tests_writer(laws,testFile)
    testFile.write(s.getvalue())
    testFile.close()

    twin_laws+=laws    
    
    
    
    
    
    
    self.reset_threshold()

    return twin_laws

  def tests_writer(self,laws,testFile):
    laws=self.purge_duplicates(laws)
    for law in laws:
      matrix=law.hkl_rotation
      basf=law.rbasf[0]
      r_red=law.rbasf[1:]
      axis=law.twin_axis
      angle=law.angle
      testFile.write(str(matrix)+"\n")
      testFile.write("Axis: "+str(axis)+"\n")
      testFile.write("Angle: "+str(angle)+"\n")
      testFile.write("Basf: "+str(basf)+"\n")
      testFile.write("R-reductions: "+str(r_red)+"\n")

  def find_twin_axes_sphere(self, hkl, threshold=None,size=10):

    model=self.model
    twin_laws=[]
    recip_orth=self.recip_orth
    metrical_matrix=self.metrical
    metrical_inverse=self.metrical_inv
    orthogonalization_inverse=self.recip_orth_inv
    a_star=recip_orth.T[0]
    b_star=recip_orth.T[1]
    c_star=recip_orth.T[2]
    b_star_cross_c_star=numpy.array(self.cross_product(b_star,c_star))
    c_star_cross_a_star=numpy.array(self.cross_product(c_star,a_star))
    a_star_cross_b_star=numpy.array(self.cross_product(a_star,b_star))
    if not threshold:
      threshold=self.cartesian_threshold

    hkl_choice=hkl[:size]

    viable_rotation_points=[]

    for i,vec0 in enumerate(hkl_choice):
      if len(twin_laws)>2*self.number_laws:
        print("Twin laws found on try " + str(i))
        twin_laws=self.purge_duplicates(twin_laws)
        if len(twin_laws)>=self.number_laws:
          break
      viable_rotation_points+=[self.find_same_distance_points(vec0,threshold)]
      viable_rotation_points_0=viable_rotation_points[i]
      v=vec0
      v_size=self.size_of_3d_vector(numpy.dot(recip_orth,v))
      #all_v_minus_w=v-viable_rotation_points_0
      v_euclid=numpy.dot(recip_orth,v)
      j=-1
      while j<i-1:
        j=j+1
        vec1=hkl_choice[j]
        if all(vec0==vec1):
          continue
        viable_rotation_points_1=viable_rotation_points[j]
        p=vec1
        p_size=self.size_of_3d_vector(numpy.dot(recip_orth,p))
        p_euclid=numpy.dot(recip_orth,p)
        for w in viable_rotation_points_0:
          w_scaled=numpy.array(w)*v_size/self.size_of_3d_vector(numpy.dot(recip_orth,w))
          for q in viable_rotation_points_1:
            q_scaled=numpy.array(q)*p_size/self.size_of_3d_vector(numpy.dot(recip_orth,q))
            cross=self.cross_product((v-w_scaled),(p-q_scaled))
            axis=cross[0]*b_star_cross_c_star+cross[1]*c_star_cross_a_star+cross[2]*a_star_cross_b_star
            if self.size_of_3d_vector(axis)<1e-16:
              continue
            unit_axis=axis/self.size_of_3d_vector(axis)
            #parallel and perpendicular components to the axis. Theoretically v_par and w_par are the same.
            v_par=numpy.dot(unit_axis,v_euclid)*unit_axis
            v_per=v_euclid-v_par
            w_per=numpy.dot(recip_orth,w_scaled)-v_par
            p_par=numpy.dot(unit_axis,p_euclid)*unit_axis
            p_per=p_euclid-p_par
            q_per=numpy.dot(recip_orth,q_scaled)-p_par

            cosvw_angle=numpy.dot(v_per,w_per)/(self.size_of_3d_vector(v_per)*self.size_of_3d_vector(w_per))
            cospq_angle=numpy.dot(p_per,q_per)/(self.size_of_3d_vector(p_per)*self.size_of_3d_vector(q_per))

            if abs(cosvw_angle)>1 or abs(cospq_angle)>1:
              continue

            vw_angle=numpy.arccos(cosvw_angle)
            pq_angle=numpy.arccos(cospq_angle)
            #take positive as anticlockwise and negative as clockwise

            if abs(vw_angle-pq_angle)<0.10472: #needs to be larger  than 0.012 so put at 1 degree (0.0175) for now - changed to 6 degrees ()
              v_w_dir=self.cross_product(v_per,w_per)
              p_q_dir=self.cross_product(p_per,q_per)
              if numpy.dot(v_w_dir,p_q_dir)<0:
                #they rotate in opposite directions
                continue
              elif numpy.dot(v_w_dir,unit_axis)<0:
                vw_angle*=-1
                pq_angle*=-1
              #the above should work but takes waaaay too long
              average_angle=(vw_angle+pq_angle)/2
              rotation=self.make_real_lattice_rotation(unit_axis, average_angle)
              rotation_lattice=numpy.dot(orthogonalization_inverse,numpy.dot(rotation,recip_orth))
              prot=numpy.dot(rotation_lattice,p)
              vrot=numpy.dot(rotation_lattice,v)
              if any(prot-q_scaled>0.1) or any(vrot-w_scaled>0.1):
                rotation_lattice=numpy.linalg.inv(rotation_lattice)
                prot=numpy.dot(rotation_lattice,p)
                vrot=numpy.dot(rotation_lattice,v)
                if any(prot-q_scaled>0.1) or any(vrot-w_scaled>0.1):
                    continue

              duplicate=False
              for twin_law in twin_laws:
                if numpy.allclose(twin_law.hkl_rotation,rotation_lattice):
                  duplicate=True
              if duplicate:
                continue

              #reciprocal axis
              rotated_hkl=numpy.dot(hkl,rotation_lattice.T)

              #if not self.sufficient_overlaps(rotated_hkl, threshold):
                #continue
              hkl_displacement=self.find_fom(rotated_hkl)

              if (hkl_displacement<threshold):
                rbasf=self.basf_estimate(rotation_lattice)
                if rbasf[0]>1e-10:
                  twin_laws+=[twin_rules("Alt",numpy.dot(orthogonalization_inverse,unit_axis),rotation_lattice,average_angle,hkl_displacement,threshold,rbasf)]


    twin_laws=self.do_rounding(twin_laws)
    twin_laws=self.purge_duplicates(twin_laws)
    return twin_laws

  def find_twofold_axes_sphere(self,hkl,threshold=None,size=10):

    model=self.model
    twin_laws=[]
    recip_orth=self.recip_orth
    metrical_matrix=self.metrical
    metrical_inverse=self.metrical_inv
    recip_orth_inv=self.recip_orth_inv
    number_laws=self.number_laws
    a_star=recip_orth.T[0]
    b_star=recip_orth.T[1]
    c_star=recip_orth.T[2]
    b_star_cross_c_star=numpy.array(self.cross_product(b_star,c_star))
    c_star_cross_a_star=numpy.array(self.cross_product(c_star,a_star))
    a_star_cross_b_star=numpy.array(self.cross_product(a_star,b_star))

    if not threshold:
      threshold=self.cartesian_threshold
    maxRdrop=0

    hkl_choice=hkl[:size]

    viable_rotation_points=[]
    #for vec in hkl_choice:
      #viable_rotation_points+=[self.find_same_distance_points(vec,threshold)]

    #R=I+sintheta K+(1-costheta )K^2
    #R=I+2K^2

    for i,v in enumerate(hkl_choice):
      if len(twin_laws)>2*number_laws:
        print("Twin laws found on try " + str(i))
        twin_laws=self.purge_duplicates(twin_laws)
        if len(twin_laws)>=number_laws:
          break
      v_size=self.size_of_3d_vector(numpy.dot(recip_orth,v))
      viable_rotation_points+=[self.find_same_distance_points(v,threshold)]
      viable_rotation_points_0=viable_rotation_points[i]
      for w in viable_rotation_points_0:
        if all(w==-v):
          continue
        w_scaled=numpy.array(w)*v_size/self.size_of_3d_vector(numpy.dot(recip_orth,w))
        axis=v+w_scaled
        axis_euclid=numpy.dot(recip_orth,axis)
        axis_unit=axis_euclid/self.size_of_3d_vector(axis_euclid)
        k=numpy.array([[0,-axis_unit[2],axis_unit[1]],[axis_unit[2],0,-axis_unit[0]],[-axis_unit[1],axis_unit[0],0]])
        k2=numpy.dot(k,k)
        rot_mat_euclid=numpy.eye(3)+2*k2
        rotation_lattice=numpy.dot(recip_orth_inv,numpy.dot(rot_mat_euclid,recip_orth))


        duplicate=False
        for twin_law in twin_laws:
          if numpy.allclose(twin_law.hkl_rotation,rotation_lattice):
            duplicate=True
        if duplicate:
          continue

        rotated_hkl=numpy.dot(hkl,rotation_lattice.T)

        #if not self.sufficient_overlaps(rotated_hkl, threshold):
          #continue
        hkl_displacement=self.find_fom(rotated_hkl)

        if (hkl_displacement<threshold):
          rbasf=self.basf_estimate(rotation_lattice)
          if rbasf[0]>1e-10 and rbasf[2]>maxRdrop/2:
            twin_laws+=[twin_rules("Alt",v+w_scaled,rotation_lattice,numpy.pi,hkl_displacement,threshold,rbasf)]
            maxRdrop=max(maxRdrop,rbasf[2])


    twin_laws=self.do_rounding(twin_laws)
    twin_laws=self.purge_duplicates(twin_laws)
    return twin_laws








  def find_same_distance_points(self, hkl, threshold=None):
    #finds points with the same d-spacing as a reciprocal lattice point hkl

    model=self.model
    rec_metrical=model.unit_cell().reciprocal_metrical_matrix()
    metrical=model.unit_cell().metrical_matrix()
    orthogonalization=self.recip_orth
    a_star=orthogonalization[0]
    b_star=orthogonalization[1]
    c_star=orthogonalization[2]
    mod_a_star=numpy.sqrt(rec_metrical[0])
    mod_b_star=numpy.sqrt(rec_metrical[1])
    mod_c_star=numpy.sqrt(rec_metrical[2])
    a_star_dot_b_star=rec_metrical[3]
    a_star_dot_c_star=rec_metrical[4]
    b_star_dot_c_star=rec_metrical[5]
    mod_a_square=metrical[0]
    mod_b_square=metrical[1]
    mod_c_square=metrical[2]
    a_dot_b=metrical[3]
    a_dot_c=metrical[4]
    b_dot_c=metrical[5]
    cos_gam_star=a_star_dot_b_star/(mod_a_star*mod_b_star)
    cos_alp_star=b_star_dot_c_star/(mod_c_star*mod_b_star)
    cos_bet_star=a_star_dot_c_star/(mod_a_star*mod_c_star)
    sin_alp_star=numpy.sqrt(1-cos_alp_star**2)

    if not threshold:
      threshold=self.cartesian_threshold    
    

    size_hkl=self.size_of_3d_vector(numpy.dot(orthogonalization,hkl))
    max_distance=size_hkl+threshold
    min_distance=size_hkl-threshold
    same_distance_points=[]

    k_coeff_ext=-mod_a_star*(cos_gam_star-cos_bet_star*cos_alp_star)/(mod_b_star*sin_alp_star**2)
    k_sqrt_mult=1/(mod_b_star*sin_alp_star)
    k_coeff_int=mod_a_star**2*(cos_gam_star**2+cos_alp_star**2+cos_bet_star**2-1-2*cos_gam_star*cos_alp_star*cos_bet_star)/sin_alp_star**2
    #The initial method to find the ranges was underestimating. Trying a new one. #now changed to squares as that's what my re-do of the maths implies!

    #h_lambda=[mod_a_square,a_dot_b,a_dot_c]
    #k_lambda=[a_dot_b,mod_b_square,b_dot_c]
    #l_lambda=[a_dot_c,b_dot_c,mod_c_square]

    #lambda_h=1/max_distance*(self.size_of_3d_vector(numpy.dot(orthogonalization,h_lambda))**2)
    #lambda_k=1/max_distance*(self.size_of_3d_vector(numpy.dot(orthogonalization,k_lambda))**2)
    #lambda_l=1/max_distance*(self.size_of_3d_vector(numpy.dot(orthogonalization,l_lambda))**2)

    #try to get the maximum h-value

    #max_h=numpy.floor(h_lambda[0]/lambda_h)
    #max_k=numpy.floor(k_lambda[1]/lambda_k)
    #max_l=numpy.floor(l_lambda[2]/lambda_l)
    
    max_h=max_distance*numpy.sqrt(mod_a_square)

    for h in range(0,int(math.floor(max_h))+1):
      k_sqrt=k_sqrt_mult*numpy.sqrt(max_distance**2+h**2*k_coeff_int)
      min_k=k_coeff_ext*h-k_sqrt
      max_k=k_coeff_ext*h+k_sqrt
      if h==0:
        min_k=max(0,min_k)
      for k in range(int(math.ceil(min_k)),int(math.floor(max_k))+1):
        discriminant_addition=-(h**2*mod_a_star**2+k**2*mod_b_star**2+2*h*k*a_star_dot_b_star)+(h*a_star_dot_c_star+k*b_star_dot_c_star)**2/mod_c_star**2
        constant_l=(h*a_star_dot_c_star+k*b_star_dot_c_star)/mod_c_star**2
        upper_discriminant=max_distance**2+discriminant_addition
        lower_discriminant=min_distance**2+discriminant_addition
        #still sometimes has upper <0 even when a twin laws is there. Needs work
        if upper_discriminant<0:
          continue
        l_U_max=-constant_l+math.sqrt(upper_discriminant)/mod_c_star#numpy.trunc(-constant_l+math.sqrt(upper_discriminant)/mod_c_star)
        l_U_min=-constant_l-math.sqrt(upper_discriminant)/mod_c_star#numpy.trunc(-constant_l-math.sqrt(upper_discriminant)/mod_c_star)
        if k==0:
          l_U_min=max(0,l_U_min)
        if lower_discriminant<0:
          possible_l=list(range(int(math.ceil(l_U_min)),int(math.floor(l_U_max)+1)))
        else:
          l_L_max=-constant_l+math.sqrt(lower_discriminant)/mod_c_star#numpy.trunc(-constant_l+math.sqrt(lower_discriminant)/mod_c_star)
          l_L_min=-constant_l-math.sqrt(lower_discriminant)/mod_c_star#numpy.trunc(-constant_l-math.sqrt(lower_discriminant)/mod_c_star)
          if k==0:
            l_L_min=max(0,l_L_min)
            if l_L_min==0:
              possible_l=list(range(int(math.ceil(l_L_max)),int(math.floor(l_U_max))+1))
            else:
              possible_l=list(range(int(math.ceil(l_U_min)),int(math.floor(l_L_min))+1))+list(range(int(math.ceil(l_L_max)),int(math.floor(l_U_max))+1))
          else:
            possible_l=list(range(int(math.ceil(l_U_min)),int(math.floor(l_L_min))+1))+list(range(int(math.ceil(l_L_max)),int(math.floor(l_U_max))+1))
        for l in possible_l:
          size_h=self.size_of_3d_vector(numpy.dot(orthogonalization,[h,k,l]))
          if hkl[0]==h and hkl[1]==k and hkl[2]==l:
            continue
          if abs(size_hkl-size_h)<threshold:
            same_distance_points+=[[h,k,l]]

    same_distance_points+=numpy.negative(same_distance_points).tolist()
    return same_distance_points

  def check_basf(self, twin_laws):
    twin_laws_2=[]
    for twin_law in twin_laws:
      if twin_law.rbasf[0]>1e-15:
        twin_laws_2+=[twin_law]
        continue
      rbasf=self.basf_estimate(twin_law.hkl_rotation)
      twin_law.rbasf=rbasf
      if rbasf[0]>1e-15:
        twin_laws_2+=[twin_law]
    return twin_laws_2

  def cross_product(self,x,y):
    z=[0,0,0]
    z[0]=x[1]*y[2]-x[2]*y[1]
    z[1]=x[2]*y[0]-x[0]*y[2]
    z[2]=x[0]*y[1]-x[1]*y[0]
    return z

  def size_of_3d_vector(self,x):
    z=0
    z+=x[0]*x[0]
    z+=x[1]*x[1]
    z+=x[2]*x[2]
    z=z**0.5
    return z

  def test_law(self,twin_law,threshold):
    new_hkl=numpy.dot(twin_law,self.bad_hkl)
    #if self.sufficient_overlaps(new_hkl,threshold)==0:
      #return None
    fom=self.find_fom()
    if fom<threshold:
      twin_law.fom=fom
    else:
      return None
    #if self.basf_estimate_short(twin_law)[0]==0:
      #return None
    basf=self.basf_estimate(twin_law)
    if basf[0]==0:
      return None
    else:
      twin_law.rbasf=basf
      return twin_law

def format_twin_string_from_law(twin_law):
  twin_str_l = []
  twin_str = "Could not get twin string"
  try:
    twin_law = twin_law.tolist()
    if len(twin_law) == 3:
      l = []
      for row in twin_law:
        for item in row:
          l.append(item)
    twin_law = l
  except:
    pass
  for ele in twin_law:
    if str(ele) == "0.0":
      ele = "0"
    elif str(ele) == "-1.0":
      ele = "-1"
    elif str(ele) == "1.0":
      ele = "1"
    else:
      ele = round(ele, 4)
    twin_str_l.append("%s" % ele)
  twin_str = " ".join(twin_str_l)
  return twin_str

def get_twinning_result_filename():
  import ntpath
  import os

  _ = ntpath.basename(OV.HKLSrc())
  _ = _.rstrip(".hkl")
  if "_twin" in _:
    _ = _.split("_twin")[0]
  _ = os.path.join(OV.StrDir(), _ + "_twinning.html")
  return _
OV.registerFunction(get_twinning_result_filename)

def on_twin_image_click(run_number):
  global twin_laws_d

  init_twin_gui()
  #if not twin_laws_d:
    #import pickle as pickle
    #p = os.path.join(OV.StrDir(), 'twin_laws_d.pickle')
    #with open(p, "rb") as infile:
      #twin_laws_d = pickle.load(infile)
    #make_law_images(twin_law_d, len(twin_laws_d))
    #for number in twin_laws_d:
      #law = twin_laws_d[number]
      #im = law.get('law_image')
      #OlexVFS.save_image_to_olex(im, "IMG_LAW%s"%number)

  try:
    twin_law = numpy.array(twin_laws_d[int(run_number)]['matrix'])
    twin_law_rnd = numpy.rint(twin_law)
    twin_law_disp_rnd = format_twin_string_from_law(twin_law_rnd)
    twin_law_disp = format_twin_string_from_law(twin_law)
    basf = float(twin_laws_d[int(run_number)]['BASF'])
    hklf5name = twin_laws_d[int(run_number)].get('hklf5name', None)
  except:
    twin_law = numpy.array(twin_laws_d[int(run_number)]['matrix'])


  if(numpy.any(numpy.abs(twin_law-twin_law_rnd)>0.05)):
    print("Using twin law: %s" %twin_law_disp)
    print("This is a non-integral twin law, and a corresponding hklf 5 format file has been made.")
    # non integral twin law, need hklf5
    OV.DelIns("TWIN")
    olx.HKLF(2)
    OV.DelIns("BASF")
    OV.AddIns("BASF %f" % basf)
    if not hklf5name:
      hklf5name = "%s_twin%02d.hkl" % (OV.HKLSrc().rsplit('\\', 1)[-1].rstrip(".hkl"), int(run_number))
    _ = os.path.join(OV.FilePath(), hklf5name)
    OV.HKLSrc(_)
  else:
    print("Using twin law: %s" %twin_law_disp)
    print("This is an integral twin law, and twinning will be handled by the refinement program.")
    OV.DelIns("TWIN")
    olx.HKLF(0)
    OV.DelIns("BASF")
    OV.AddIns("BASF %f"%basf)
    OV.AddIns("TWIN %s" % twin_law_disp_rnd)
  OV.UpdateHtml()
OV.registerFunction(on_twin_image_click)

def write_twin_images_to_disk(name, fn_base):
  l = ['on', 'off']
  for state in l:
    n = '%s%s.png'%(name,state)
    im_data = OlexVFS.read_from_olex(n)
    image_name = fn_base+n
    with open (image_name, 'wb') as wFile:
      wFile.write(im_data)
  return image_name

def reset_twin_law_img():
  global twin_laws_d
  olex_refinement_model = OV.GetRefinementModel(False)
  if 'twin' in olex_refinement_model:
    c = olex_refinement_model['twin']['matrix']
    curr_law = []
    for row in c:
      for el in row:
        curr_law.append(el)
    for i in range(3):
      curr_law.append(0.0)
    curr_law = tuple(curr_law)

  else:
    curr_law = (1, 0, 0, 0, 1, 0, 0, 0, 1)
  for law in twin_laws_d:
    name = twin_laws_d[law]['name']
    matrix = twin_laws_d[law]['matrix']
    if curr_law == matrix:
      OV.CopyVFSFile("%son.png" %name, "%s.png" %name,2)
    else:
      OV.CopyVFSFile("%soff.png" %name, "%s.png" %name,2)
  OV.UpdateHtml()
OV.registerFunction(reset_twin_law_img)

OlexCctbxTwinLaws_instance = OlexCctbxTwinLaws()

#taken from stackoverflow 8391411
class HiddenPrints:
  def __enter__(self):
    self._original_stdout = sys.stdout
    sys.stdout = open(os.devnull, 'w')

  def __exit__(self, exc_type, exc_val, exc_tb):
    try:
      sys.stdout.close()
      sys.stdout = self._original_stdout
    except:
      pass


def round_matrix_elements(matrix):
  round_tol = 0.01
  round_val_l = [0,0.25,0.5,0.75,1,0.333,0.666]
  round_disp_l = ['0','1/4','1/2','3/4','1','1/3','2/3']
  round_l = list(zip(round_val_l, round_disp_l))
  new_matrix = []
  for item in matrix:
    try:
      item = round(item,3)
    except:
      m.append(item)
      continue
    sign = ""
    if abs(item) != item:
      sign = "-"
    if str(item) == "0.0":
      item = "0"
    elif str(item) == "1.0":
      item = "1"
    elif str(item) == "-1.0":
      item = "-1"
    else:
      for val,disp in round_l:
        if val-round_tol < abs(item) < val+round_tol:
          item = sign+disp
          break
    new_matrix.append(item)
  return new_matrix



def find_2_fold_olex2():
  import gui.tools
  cmd = "TestR"
  res = [_f for _f in gui.tools.scrub(cmd) if _f]

  if "matrix" not in " ".join(res).lower():
    return
  matrix = res[-3:]
  i = 0
  for line in res:
    if line.startswith('wR2'):
      basf = res[i+1].split()[2].strip()
    i += 1
  basf = 0.2
  r = 0.05
  r_diff = -0.03
  d = {}
  d.setdefault(1, {})
  matrix = [float(i) for i in " ".join(matrix).split()]
  d[1]['matrix'] = matrix
  d[1]['BASF'] = basf
  d[1]['r'] = r
  d[1]['r_diff'] = r_diff
  d[1]['name'] = "law%s" % 1
  d[1]['HKLSrc'] = OV.HKLSrc()
  d[1]['hklf5name'] = "%s.olex2_hklf5.hkl" % OV.FileName()

  make_twinning_gui(d)
OV.registerFunction(find_2_fold_olex2)

def make_law_images(twin_law, lawcount):
  global twin_laws_d
  from ImageTools import IT
#  IT = ImageTools()
  from PilTools import MatrixMaker
  MM = MatrixMaker()

  bg_col = OV.GetParam('gui.html.table_bg_colour')
  bg_col_img = IT.RGBToHTMLColor(IT.adjust_colour(bg_col.rgb,luminosity=1.3))
  font_col = "#444444"
  font_col_basf = "#447744"

  basf = twin_law['BASF']
  r = twin_law['r']
  r_diff = twin_law['r_diff']
  matrix = twin_law['matrix']
  name = "law%s" %lawcount
  if basf == "n/a":
    font_col_basf = OV.GetParam('gui.blue').rgb
  elif float(basf) < 0.1:
    font_col_basf = OV.GetParam('gui.red').rgb
    basf = "%.2f" %float(basf)
  else:
    font_col_basf = OV.GetParam('gui.green_text').rgb
    basf = "%.2f" %float(basf)

  txt = [{'txt':"R=%.2f%%, -%.2f%%" %((float(r)*100),(float(r_diff)*100)),
          'font_colour':font_col},
         {'txt':"basf=%s" %str(basf),
          'font_colour':font_col_basf}]
  states = ['on', 'off', 'hover', 'down']
  rounded_matrix = round_matrix_elements(matrix)
  for state in states:
    image_name, img  = MM.make_3x3_matrix_image(name, rounded_matrix, txt, state, bar_col=font_col_basf, bgcolor=bg_col_img)
    twin_laws_d[lawcount].setdefault('law_image_name', image_name)

def init_twin_gui():
  global twin_laws_d
  if twin_laws_d:
    if not twin_laws_d[1]['HKLSrc'] == OV.HKLSrc():
      twin_laws_d = {}

  if not twin_laws_d:
    import pickle as pickle
    p = os.path.join(OV.StrDir(), 'twin_laws_d.pickle')
    if os.path.exists(p):
      with open(p, "rb") as infile:
        twin_laws_d = pickle.load(infile)
    else:
      twin_laws_d = {}

  if twin_laws_d:
    for idx in twin_laws_d:
      make_law_images(twin_laws_d[idx], len(twin_laws_d))
    #for number in twin_laws_d:
      #law = twin_laws_d[number]
      #im = law.get('law_image')
      #OlexVFS.save_image_to_olex(im, "IMG_LAW%s" % number)
OV.registerFunction(init_twin_gui, False, 'twin')

def make_twinning_gui(laws):
  font_col = "#444444"
  font_col_basf = "#447744"
  bg_col = OV.GetParam('gui.html.table_bg_colour')

  global twin_laws_d
  twin_laws_d = laws
  p = os.path.join(OV.StrDir(), 'twin_laws_d.pickle')
  import pickle as pickle
  pickle.dump(twin_laws_d, open(p, "wb"))

  r_list = []
  for count in laws:
    twin_law_d = laws[count]
    make_law_images(twin_law_d, count)
    r_list.append((twin_law_d['r'], count, twin_law_d['BASF']))

  history = ""
  ins_file = ""
  law_txt = ""
  i = 0
  html = "<tr bgcolor='%s'><td></td><td>" %bg_col
  fn_base = get_twinning_result_filename().rstrip(".html")

  for r, run, basf in r_list:
    i += 1
    image_name = twin_laws_d[run].get('law_image_name', "XX")
    #write_twin_images_to_disk(image_name, fn_base)
    use_image = "%s%soff.png" %(fn_base, image_name)
    img_src = "%s.png" %image_name
    name = twin_laws_d[run].get('name', "XX")
    #href = 'spy.on_twin_image_click(%s)'
    href = 'spy.on_twin_image_click(%s)>>html.Update' % (i,)
    law_txt = "<a href='%s'><img src=%s></a>&nbsp;" %(href, use_image)

    d = {'name':image_name,
         'nameupper':image_name.upper(),
         'tool_img':image_name,
         'down':'down',
         'up':'up',
         'on':'on',
         'hover':'hover',
         'cmds':href,
         'target':"",
         'feedback':"",
         'bgcolor':"#ff9999",
    }

    law_txt = '''
    <font size='$GetVar(HtmlFontSizeControls)'>
    <input
    name="IMG_%(nameupper)s"
    type="button"
    image="up=%(tool_img)s%(on)s.png,down=%(tool_img)s%(down)s.png,hover=%(tool_img)s%(hover)s.png"
    hint="%(target)s"
    onclick="%(cmds)s%(feedback)s"
    bgcolor="%(bgcolor)s"
  >
  </font>
  '''%d
    #self.twin_law_gui_txt += "%s" %(law_txt)
    control = "IMG_%s" %image_name.upper()
    html += law_txt

  #twin_laws_d[lawcount] = {'number':lawcount,
                                #'law':matrix,
                                #'R1':r,
                                #'BASF':basf,
                                #'law_image':img,
                                #'law_txt':law_txt,
                                #'law_image_name':image_name,
                                #'name':name,
                                #'ins_file':f_data,
                                #'history':history,
                                #}
  #l += 1


  html += "</td></tr>"
  write_twinning_result_to_disk(html)
  OV.UpdateHtml()

def write_twinning_result_to_disk(html):
  with open(get_twinning_result_filename(), 'w') as wFile:
    wFile.write(html)
