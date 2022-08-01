
from olexFunctions import OlexFunctions
OV = OlexFunctions()

import os
import shutil
import htmlTools
import olex
import olx
import gui
import numpy as np
import sys
import copy
import scipy as sp

import boost.python
ext = boost.python.import_ext("smtbx_refinement_least_squares_ext")

import time
debug = bool(OV.GetParam("olex2.debug", False))

import NoSpherA2
import shutil
import leverage
from refinement import FullMatrixRefine

from scitbx.array_family import flex
from cctbx import xray

from cctbx import adptbx

from leverage import leverage_normal_eqns

#n=NoSpherA2()
#n.launch() #NoSpherA2.NoSpherA2.NoSpherA2_instance.launch()

from cctbx_olex_adapter import OlexCctbxAdapter

instance_path = OV.DataDir()


np.set_printoptions(threshold=sys.maxsize)



try:
  from_outside = False
  p_path = os.path.dirname(os.path.abspath(__file__))
except:
  from_outside = True
  p_path = os.path.dirname(os.path.abspath("__file__"))

l = open(os.sep.join([p_path, 'def.txt'])).readlines()
d = {}
for line in l:
  line = line.strip()
  if not line or line.startswith("#"):
    continue
  d[line.split("=")[0].strip()] = line.split("=")[1].strip()

p_name = d['p_name']
p_htm = d['p_htm']
p_img = eval(d['p_img'])
p_scope = d['p_scope']

OV.SetVar('LauraNonSpherTests_plugin_path', p_path)

from PluginTools import PluginTools as PT

#refine_extinction=True #True #False

class LauraNonSpherTests(PT):

  def __init__(self):
    super(LauraNonSpherTests, self).__init__()
    self.p_name = p_name
    self.p_path = p_path
    self.p_scope = p_scope
    self.p_htm = p_htm
    self.p_img = p_img
    self.deal_with_phil(operation='read')
    self.print_version_date()
    if not from_outside:
      self.setup_gui()
    OV.registerFunction(self.make_gradient_file,True,"LauraNonSpherTests")
    OV.registerFunction(self.rename_gradient_file,True,"LauraNonSpherTests")
    OV.registerFunction(self.rename_comparison_file,True,"LauraNonSpherTests")
    OV.registerFunction(self.do_test,True,"LauraNonSpherTests")
    OV.registerFunction(self.refine,True,"LauraNonSpherTests")
    OV.registerFunction(self.refine_multiple,True,"LauraNonSpherTests")
    OV.registerFunction(self.epsilon_tests,True,"LauraNonSpherTests")
    OV.registerFunction(self.start_data,True,"LauraNonSpherTests")
    OV.registerFunction(self.deposit_phrase,True,"LauraNonSpherTests")
    OV.registerFunction(self.spherical_refinement_with_ns_fcalc,True,"LauraNonSpherTests")
    OV.registerFunction(self.make_hkl_from_data,True,"LauraNonSpherTests")
    OV.registerFunction(self.make_fc_hkl,True,"LauraNonSpherTests")

    # END Generated =======================================

  def make_gradient_file(self):
    #creates a 'nospte' instance and asks it to make the gradient file, which is labelled 'name_gradient_records.txt', and contains, in order: 
    #D_num
    #D_appr(num)
    #D_appr(ana)
    #D_hyb
    #D_sph(num)
    #D_sph(ana)
    #Fc^2_appr
    #W_appr
    #Fc^2_sph
    #W_sph
    #Orthogonalization Matrix 'A'
    #Current atomic parameters
    #Weights retrieved differently

    OV.SetParam('snum.NoSpherA2.no_backup',True)
    nospte = NonSphereTests()
    nospte.run()
    OV.SetParam('snum.NoSpherA2.no_backup',False)
    nospte = None

  def rename_gradient_file(self,append):
    #utilised by testing functions to rename gradient files for storage
        #    shutil.copyfile((os.path.split(OV.HKLSrc())[-1]).rsplit('.',1)[0]+'.tsc', 'originaltsc.tsc')
    filename=(os.path.split(OV.HKLSrc())[-1]).rsplit('.',1)[0]+'_gradient_records.txt'    
    os.rename(filename, (os.path.split(OV.HKLSrc())[-1]).rsplit('.',1)[0]+'_gradient_records_'+append+'.txt' )
  def rename_comparison_file(self,name=""):
    #utilised by testing functions to rename the comparisonData file to store specific details in specifically names files
    if name=="":
      newname="comparisonData_"+(os.path.split(OV.HKLSrc())[-1]).rsplit('.',1)[0]+".txt"
    else:
      newname=name+"_"+(os.path.split(OV.HKLSrc())[-1]).rsplit('.',1)[0]+".txt"
    os.rename("comparisonData.txt",newname)

  def refine(self):
    #performs one numerical refinement step
    nospte=NonSphereTests()
    nospte.refinement_step()
    nospte = None

  def refine_multiple(self,hybrid=False):
    #performs 20 numerical refinement steps and records them in the comparisonData file. At the end it returns the structure which obtained the lowest wR2 value
    OV.SetParam('snum.NoSpherA2.no_backup',True)
    nospte=NonSphereTests()
    nospte.run(keep_old_files=False)
    nospte = None
    bestStructure=storedStructure(None, 1.)
    i=0
    for i in range(5):
      nospte=NonSphereTests()
      nospte.run(keep_old_files=True)
      bestStructure=nospte.compare_structures(bestStructure)[0]
      nospte.get_full_data(steps=i,hybrid=hybrid)
      nospte.refinement_step(hybrid=hybrid)
      nospte = None
    nospte = None
    nospte=NonSphereTests()
    nospte.run(keep_old_files=True)
    bestStructure,currentIsBest=nospte.compare_structures(bestStructure)
    nospte.get_full_data(steps=20,hybrid=hybrid)
    if not currentIsBest:
      nospte.replace_structure(bestStructure)
    nospte.run(keep_old_files=True)
    nospte.deposit_phrase("Final Structure")
    nospte.get_full_data(steps="Final",hybrid=hybrid)
    nospte = None
    shutil.rmtree(r"olex2\Wfn_job")
    shutil.rmtree(r"olex2\NoSpherA2_history")

    OV.SetParam('snum.NoSpherA2.no_backup',False)

  def spherical_refinement_with_ns_fcalc(self):
    #performs 'mixed' refinement, which uses spherical derivatives but nonspherical Fcalcs
    #this performs almost as well as approximate nonspherical refinement, but takes the same time, so we do not see any practical applications
    #like 'refine_multiple', this performs 20 steps and stores them in 'comparisonData'
    OV.SetParam('snum.NoSpherA2.no_backup',True)
    OV.SetParam('snum.NoSpherA2.use_aspherical',False)
    olex.m('anis -h')
    for i in range(100):
      olex.m('refine')

    bestStructure=storedStructure(None, 1.)

    #put NoSpherA2 settings
    OV.SetParam('snum.NoSpherA2.use_aspherical',True)
    OV.SetParam('snum.NoSpherA2.full_HAR',False)  
    #OV.SetParam('snum.NoSpherA2.source',"ORCA") #ZZZZ Removed set ORCA
    OV.SetParam('snum.NoSpherA2.Calculate',True)  
    #OV.SetParam('snum.NoSpherA2.basis_name',"def2-TZVPP")
    OV.SetParam('snum.NoSpherA2.becke_accuracy',"High")
    OV.SetParam('snum.NoSpherA2.method',"PBE")    


    nospte=NonSphereTests()
    nospte.run()    
    nospte = None

    for i in range(20):

      nospte=NonSphereTests()
      NoSpherA2.NoSpherA2.NoSpherA2_instance.launch()
      bestStructure=nospte.compare_structures(bestStructure)[0]
      nospte.get_full_data(steps=i,mix=True)
      nospte.sph_ns_refine()  
      nospte = None


    nospte = None
    nospte=NonSphereTests()
    nospte.run(keep_old_files=True)
    bestStructure,currentIsBest=nospte.compare_structures(bestStructure)
    nospte.get_full_data(steps=20,mix=True)
    if not currentIsBest:
      nospte.replace_structure(bestStructure)
    nospte.run(keep_old_files=True)
    nospte.deposit_phrase("Final Structure")
    nospte.get_full_data(steps="Final",mix=True)
    nospte = None
    OV.SetParam('snum.NoSpherA2.no_backup',False)


  def start_data(self):
    #deposits some basic descriptive data of the tests being run
    OV.SetParam('snum.NoSpherA2.no_backup',True)
    nospte=NonSphereTests()
    nospte.run(keep_old_files=False)  
    nospte.get_base_data()
    nospte=None    
    OV.SetParam('snum.NoSpherA2.no_backup',False)

  def epsilon_tests(self):
    #generates gradient records for epsilon from 1e-1 to 1e-10. This isn't really substantially useful, mostly eliminiting those at the extreme range (1e-7 plus)
    OV.SetParam('snum.NoSpherA2.no_backup',True)
    OV.SetParam('snum.NoSpherA2.use_aspherical',False)
    olex.m('anis -h')
    for i in range(30):
      olex.m('refine')
    print("Done sph refine")
    #put NoSpherA2 settings
    OV.SetParam('snum.NoSpherA2.use_aspherical',True)
    OV.SetParam('snum.NoSpherA2.full_HAR',False)  
    #OV.SetParam('snum.NoSpherA2.source',"ORCA") #ZZZZ Removed set ORCA
    OV.SetParam('snum.NoSpherA2.Calculate',True)  
    #OV.SetParam('snum.NoSpherA2.basis_name',"def2-TZVPP")
    OV.SetParam('snum.NoSpherA2.becke_accuracy',"High")
    OV.SetParam('snum.NoSpherA2.method',"PBE")
    print("Done settings")
    olex.m('spy.LauraNonSpherTests.make_gradient_file()')   
    print("Done basic gradient file")

    import os
    olex.m('spy.LauraNonSpherTests.make_gradient_file()')   
    print("Done basic gradient no.2")
    epsilons=[1e-3,1e-4,1e-2,1e-1,1e-5,1e-6,1e-7,1e-8,1e-9,1e-10]#[1e-10,1e-9,1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1]#[1,1e-1,1e-2,5e-3,2e-3,1e-3,5e-4,2e-4,1e-4,5e-5,2e-5,1e-5,1e-6,1e-7,1e-8,1e-9,1e-10]#[1e-1,1e-2,1e-3,1e-4,1e-5,1e-6,1e-7,1e-8,1e-9,1e-10]
    print("Done eps set")
    if True:
      for eps in epsilons:
        print("doing eps"+str(eps))
        OV.SetParam('lauranonsphertests.vars.epsilon', eps)
        olex.m('spy.LauraNonSpherTests.make_gradient_file()')   
        #filename=(os.path.split(OV.HKLSrc())[-1]).rsplit('.',1)[0]+'_gradient_records.txt'    
        parameter='spy.LauraNonSpherTests.rename_gradient_file("eps_+'+str(eps)+'")'
        olex.m(parameter)
        print("done eps"+str(eps))
        #os.rename(filename, (os.path.split(OV.HKLSrc())[-1]).rsplit('.',1)[0]+'_gradient_records_eps_'+str(eps)+'.txt' ) #fv_ five
      print("Finished sph eps tests")
    OV.SetParam('lauranonsphertests.vars.epsilon', 1e-4)
    OV.SetParam('snum.NoSpherA2.no_backup',False)
    olex.m('anis -h')
    if True:
      for i in range(10):
        olex.m('refine')
        OV.SetParam('snum.NoSpherA2.Calculate',True)      
      print("Done app refine")
      for eps in epsilons:
        print("doing eps"+str(eps))
        OV.SetParam('lauranonsphertests.vars.epsilon', eps)
        olex.m('spy.LauraNonSpherTests.make_gradient_file()')   
        #filename=(os.path.split(OV.HKLSrc())[-1]).rsplit('.',1)[0]+'_gradient_records.txt'    
        parameter='spy.LauraNonSpherTests.rename_gradient_file("app_eps_+'+str(eps)+'")'
        olex.m(parameter)
        #olex.m('spy.LauraNonSpherTests.rename_gradient_file("app_eps_"+'+str(eps))
        print("done eps"+str(eps))
        #os.rename(filename, (os.path.split(OV.HKLSrc())[-1]).rsplit('.',1)[0]+'app_gradient_records_eps_'+str(eps)+'.txt' ) #fv_ five
    OV.SetParam('lauranonsphertests.vars.epsilon', 1e-3)
    OV.SetParam('snum.NoSpherA2.no_backup',False)


  def do_test(self):
    #allows simple retrieval of data like positions and such
    #ensure you 'run' the gradient records before calling this
    nospte=NonSphereTests()
    nospte.get_full_data()
    nospte = None

  def deposit_phrase(self, phrase):
    #puts a phrase into the comparisonData file
    nospte=NonSphereTests()
    nospte.deposit_phrase(phrase)
    nospte = None
  def make_hkl_from_data(self):
    #prints a hkl-ish formatted list in the same order as stored in olex2, which may not match that in the associated .hkl file.
    nospte=NonSphereTests()
    nospte.make_own_hkl()
    nospte = None
    
  def make_fc_hkl(self):
    nospte=NonSphereTests()
    nospte.make_own_hkl(fcalcs=True)
    nospte = None
    
  
LauraNonSpherTests_instance = LauraNonSpherTests()
print("OK.")

class storedStructure():
  #a small class to store the relevant details of a structure along the refinement path
  def __init__(self,atoms,wR2):
    self.atoms=atoms
    self.wR2=wR2

class NonSphereTests(OlexCctbxAdapter):
  def __init__(self):
    #sets up a lot of details such as the metrical matrix
    super(NonSphereTests, self).__init__()
    self.model=self.xray_structure() 
    self.orthogonalization=np.reshape(np.array(self.model.unit_cell().orthogonalization_matrix()),(3,3))
    self.fmr=FullMatrixRefine()
    initial_fmr=self.fmr.run(build_only=True)
    self.fmr.reparametrisation.linearise()
    reparam=self.fmr.reparametrisation
    jt=reparam.jacobian_transpose_matching_grad_fc()
    self.jacobian_transpose_matching_grad_fc=jt.as_dense_matrix().as_numpy_array()
    self.fivePointDerivative=False ##five ZZZ
    self.refine_extinction=OV.GetParam('lauranonsphertests.vars.refineextinction')
    self.grow=OV.GetParam('lauranonsphertests.vars.grow_molecule')
    self.grow=False #override
    gl=self.model.unit_cell().metrical_matrix()
    metrical = np.zeros((3,3))
    for i in range(3):
      metrical[i,i] = gl[i]
    metrical[0,1] = gl[3]
    metrical[0,2] = gl[4]
    metrical[1,2] = gl[5]
    metrical[1,0] = gl[3]
    metrical[2,0] = gl[4]
    metrical[2,1] = gl[5]    
    self.metrical=metrical
    self.depositFile="comparisonData.txt"

  def run(self,keep_old_files=False,hybrid=False):
    #generates the gradient file for the current structure
    #This is labelled 'name_gradient_records.txt', and contains, in order: 
    #D_num
    #D_appr(num)
    #D_appr(ana)
    #D_hyb
    #D_sph(num)
    #D_sph(ana)
    #Fc^2_appr
    #W_appr
    #Fc^2_sph
    #W_sph
    #Orthogonalization Matrix 'A'
    #Current atomic parameters
    #Weights retrieved differently    
    #
    #
    #This is obtained by using the jacobean to implement the movement of each parameter of the structure x in turn, and using that to calculate the derivative


    keep_old_files=False #ZZZZ Override
    epsilon=OV.GetParam('lauranonsphertests.vars.epsilon')
    model = self.model  
    fmr=FullMatrixRefine()
    initial_fmr=fmr.run(build_only=True)
    fmr.reparametrisation.linearise()
    jt=self.jacobian_transpose_matching_grad_fc

    filename = olx.FileName()


    orthogonalization=np.reshape(np.array(model.unit_cell().orthogonalization_matrix()),(3,3))
    orthogonalization_inv=np.linalg.inv(orthogonalization)  
    metrical=model.unit_cell().metrical_matrix()
    unit_cell_lengths=np.sqrt(metrical[:3])

    #this will hide the .hkl file in Wfn_job, so NoSpherA2 will not make a backup dir and tehrefore force ORCA to restart using the old wavefunction
    if keep_old_files == True:
      path = os.path.join(OV.FilePath(),"olex2","Wfn_job")
      try:
        shutil.move(os.path.join(path,olx.FileName()+".hkl"),os.path.join(path,olx.FileName()+".hidden_hkl"))
      except:
        pass    
    NoSpherA2.NoSpherA2.NoSpherA2_instance.launch()
    shutil.copyfile((os.path.split(OV.HKLSrc())[-1]).rsplit('.',1)[0]+'.tsc', 'originaltsc.tsc')
    #shutil.copyfile((OV.HKLSrc().rsplit('\\',1)[-1]).rsplit('.',1)[0]+'.tsc', 'originaltsc.tsc')
    #os.rename((OV.HKLSrc().rsplit('\\',1)[-1]).rsplit('.',1)[0]+'.tsc','originaltsc.tsc')
    tsc_design=get_design()
    fcalc_tru=self.get_fc_abs()**2
    fcalc_spher=self.get_fc_abs(table=False)**2

    f_sph=self.get_fc_abs(table=False)        

    atoms = model.scatterers()
    atom_count=len(atoms)

    shift_x=np.array([epsilon,0,0])/unit_cell_lengths
    shift_y=np.array([0,epsilon,0])/unit_cell_lengths
    shift_z=np.array([0,0,epsilon])/unit_cell_lengths

    initial_design=get_design(table=False)
    index_count=(np.array(initial_design.f_calc())).size
    parameter_count=(np.array(initial_design.design_matrix())).size//index_count

    design_ma=np.array(initial_design.design_matrix()).reshape(index_count,parameter_count)

    tsc_design_ma=np.array(tsc_design.design_matrix()).reshape(index_count,parameter_count)
    labels=parameter_labels(fmr, parameter_count)
    if not self.refine_extinction and labels[-1]=='EXTI':
      design_ma=design_ma[:,:-1]
      tsc_design_ma=tsc_design_ma[:,:-1]  
      parameter_count=parameter_count-1
    if labels[-1]=='EXTI':
      print("Extinction not yet implemented!")
      return

    y_parameter_count=0
    y_params_per_atom=[]
    for atom in atoms:
      atom_count=0
      atom_count+=3
      if atom.flags.use_u_iso():
        atom_count+=1
      elif atom.flags.use_u_aniso():
        atom_count+=6
      y_parameter_count+=atom_count
      y_params_per_atom+=[atom_count]


    grads_spher=np.full((index_count,y_parameter_count), np.nan)
    grads_aspher_tru=np.full((index_count,y_parameter_count), np.nan)
    grads_aspher_app=np.full((index_count,y_parameter_count), np.nan)
    grads_hybrid=np.full((index_count,y_parameter_count), np.nan)    


    grads_tru_second=np.full((index_count,parameter_count), np.nan)
    grads_app_second=np.full((index_count,parameter_count), np.nan)
    grads_hyb_second=np.full((index_count,parameter_count), np.nan)
    grads_sph_second=np.full((index_count,parameter_count), np.nan)
    for j,line in enumerate(jt):
      shift_indices=np.where(line==1)[0]
      atomic_shift_indices=copy.deepcopy(shift_indices)
      k=shift_indices[0]
      relevantAtom=None
      for i in range(len(atoms)):
        if k<y_params_per_atom[i]:
          relevantAtom=i
          break
        else:
          k-=y_params_per_atom[i]
          atomic_shift_indices-=y_params_per_atom[i]
      if shift_indices[-1]>np.sum(y_params_per_atom[:i+1]):
        print("Problem:"
              "Mutiple atoms per x parameter. Cannot currently evaluate")
        return
      if atomic_shift_indices[0]<3:
        #position shift
        direction=[0.,0.,0.]
        dimensions=0
        if 0 in atomic_shift_indices:
          direction+=shift_x
          dimensions+=1
        if 1 in atomic_shift_indices:
          direction+=shift_y
          dimensions+=1
        if 2 in atomic_shift_indices:
          direction+=shift_z
          dimensions+=1
        if dimensions==2:
          direction/=np.sqrt(2)
        elif dimensions==3:
          direction/=np.sqrt(3)
        gt,ga,gs=self.get_derivative(relevantAtom,direction,keep_old_files=keep_old_files,hybrid=hybrid)
        grads_tru_second[:,j]=gt
        grads_app_second[:,j]=ga
        if atoms[relevantAtom].label[0]=='H':
          grads_hyb_second[:,j]=gt
        else:
          grads_hyb_second[:,j]=ga
        grads_sph_second[:,j]=gs
        if np.any(atomic_shift_indices-2>0):
          print("Problem: Positions and adps per x parameter. Cannot currently evaluate")
          return      
      else:
        #adp shift only
        grads_tru_second[:,j]=tsc_design_ma[:,j]
        grads_app_second[:,j]=tsc_design_ma[:,j]
        grads_hyb_second[:,j]=tsc_design_ma[:,j]
        grads_sph_second[:,j]=design_ma[:,j]




    tsc_design_ma_extended=np.dot(tsc_design_ma,jt)
    design_ma_extended=np.dot(design_ma,jt)

    col=0

    grads_aspher_app=grads_app_second
    grads_aspher_tru=grads_tru_second
    grads_spher=grads_sph_second
    grads_hybrid=grads_hyb_second




    filename=(os.path.split(OV.HKLSrc())[-1]).rsplit('.',1)[0]+'_gradient_records.txt'

    record=open(filename,'w')

    record.write("Gradient Records")

    labels=parameter_labels(fmr, parameter_count)
    record.write(str(labels))
    record.write("\n")    

    if labels[-1]=='EXTI' and self.refine_extinction:
      grads_aspher_app[:,col]=tsc_design_ma[:,col]
      grads_aspher_tru[:,col]=tsc_design_ma[:,col]
      grads_hybrid[:,col]=tsc_design_ma[:,col]
      grads_spher[:,col]=design_ma[:,col]
      exti_value=fmr.run(build_only=True, normal_equations_class=leverage.leverage_normal_eqns).reparametrisation.extinction.value
      col+=1

    record.write("True Gradient with tsc recalculation \n")
    for line in grads_aspher_tru:
      record.write(np.array2string(line,max_line_width=1e10,separator=', ',formatter={'float_kind':lambda x: "%.20f" % x}))
      record.write("\n")
    record.write("\n") 

    record.write("Approximate Gradient with original tsc \n")
    for line in grads_aspher_app:
      record.write(np.array2string(line,max_line_width=1e10,separator=', ',formatter={'float_kind':lambda x: "%.20f" % x}))
      record.write("\n")
    record.write("\n")

    record.write("Tsc-based design matrix\n") 
    for line in tsc_design_ma:
      record.write(np.array2string(line,max_line_width=1e10,separator=', ',formatter={'float_kind':lambda x: "%.20f" % x}))
      record.write("\n")
    record.write("\n")    

    record.write("Hybrid Gradient with only hydrogens true \n")
    for line in grads_hybrid:
      record.write(np.array2string(line,max_line_width=1e10,separator=', ',formatter={'float_kind':lambda x: "%.20f" % x}))
      record.write("\n")
    record.write("\n")      

    record.write("Spherical Calculated Gradient\n")
    for line in grads_spher:
      record.write(np.array2string(line,max_line_width=1e10,separator=', ',formatter={'float_kind':lambda x: "%.20f" % x}))
      record.write("\n")
    record.write("\n")     

    record.write("Spherical Design\n")
    for line in design_ma:
      record.write(np.array2string(line,max_line_width=1e10,separator=', ',formatter={'float_kind':lambda x: "%.20f" % x}))
      record.write("\n")
    record.write("\n")     

    recordname,hklFile=self.recordname()
    weights=np.array(initial_design.weights())

    record.write("Fcalc true\n") #fcalc_tru
    record.write(np.array2string(fcalc_tru,max_line_width=1e10,separator=', ',formatter={'float_kind':lambda x: "%.20f" % x}))
    record.write("\n")
    record.write("\n")
    record.write("Weight true\n")
    record.write(np.array2string(weights,max_line_width=1e10,separator=', ',threshold=1000000,formatter={'float_kind':lambda x: "%.20f" % x}))
    record.write("\n")
    record.write("\n")
    record.write("Fcalc Spherical\n") #fcalc_spher
    record.write(np.array2string(fcalc_spher,max_line_width=1e10,separator=', ',formatter={'float_kind':lambda x: "%.20f" % x}))
    record.write("\n")
    record.write("\n")    
    record.write("Weight Spherical\n")
    record.write(np.array2string(weights,max_line_width=1e10,separator=', ',threshold=1000000,formatter={'float_kind':lambda x: "%.20f" % x}))

    record.write("\n")
    record.write("Orthogonalization Matrix\n")
    record.write(str(orthogonalization))
    record.write("\n")

    record.write("Parameters\n")
    model = self.model  
    atoms = model.scatterers()  
    for atom in atoms:
      atomParams=[atom.label,atom.site]
      if atom.flags.use_u_iso():
        atomParams+=[atom.u_iso]
      elif atom.flags.use_u_aniso():
        atomParams+=[atom.u_star]
      atomParams+=["Euclidean Positions:", np.dot(self.orthogonalization,atom.site)]
      record.write(str(atomParams))
      record.write("\n") 
    if labels[-1]=='EXTI' and self.refine_extinction:
      record.write("EXTI Parameter"+str(exti_value)+"\n")



    record.write("Weights Old Spher\n")
    record.write(np.array2string(np.array(initial_design.weights()),max_line_width=1e10,separator=', ',threshold=1000000,formatter={'float_kind':lambda x: "%.10f" % x}))
    record.write("\n")    

    record.write("Weights Old True\n")
    record.write(np.array2string(np.array(tsc_design.weights()),max_line_width=1e10,separator=', ',threshold=1000000,formatter={'float_kind':lambda x: "%.10f" % x}))
    record.write("\n")        
    record.close()    


    NoSpherA2.NoSpherA2.NoSpherA2_instance.launch()    

  def get_derivative(self,atomno,direction,keep_old_files=False,hybrid=False):
    #by moving in 'plus direction' and 'minus direction', we generate an approximation to the derivative through a function similar to (f(x+h)-f(x-h))/2h
    #there is still a bodge on the direction where we assume that if more than one parameter is present it is because the direction is of the style (1, 1, 0) and that these will never be
    #not the same size. The 'unit' vectors as per the jacobean could be [1,0,0] or [1,1,0], either of which are divided by a size of '1', the second just probably has a greater impact.


    keep_old_files=False #ZZZ OVERRIDE
    hybrid_skip=False
    #direction in relative coordinates please!
    fullsize=np.max(direction)

    atoms = self.model.scatterers()
    atom=atoms[atomno]
    atomtype=atom.element_symbol()
    if hybrid and atomtype!='H':
      hybrid_skip=True
    old_location=atom.site   

    new_location_plus=old_location+direction
    new_location_minus=old_location-direction


    #this will hide the .hkl file in Wfn_job, so NoSpherA2 will not make a backup dir and tehrefore force ORCA to restart using the old wavefunction
    if keep_old_files == True:
      path = os.path.join(OV.FilePath(),"olex2","Wfn_job")
      try:
        shutil.move(os.path.join(path,olx.FileName()+".hkl"),os.path.join(path,olx.FileName()+".hidden_hkl"))
      except:
        pass         

    fcx0=self.get_fc_complex(original=True)   
    fcxsph_0=self.get_fc_complex(table=False)      



    id = self.olx_atoms.atom_ids[atomno]

    if self.grow:
      olex.m('fuse')

    atom.site=new_location_plus
    olx.xf.au.SetAtomCrd(id, *new_location_plus)    # olx.xf.au.SetAtomU(id, *u_trans)
    olx.xf.EndUpdate()    
    #First fcalc

    if self.grow:
      olex.m('pack cell')
      olex.m('grow')

    f1_plus = self.get_fc_abs(original=True)  
    fcx1_plus=self.get_fc_complex(original=True)
    f_spher_plus= self.get_fc_abs(table=False)  
    fcxsph_plus=self.get_fc_complex(table=False)


    if self.grow:
      olex.m('fuse')

    atom.site=new_location_minus
    olx.xf.au.SetAtomCrd(id, *new_location_minus)    
    olx.xf.EndUpdate() 

    if self.grow:
      olex.m('pack cell')
      olex.m('grow')

    f1_minus= self.get_fc_abs(original=True)  
    fcx1_minus=self.get_fc_complex(original=True)
    f_spher_minus= self.get_fc_abs(table=False)    
    fcxsph_minus=self.get_fc_complex(table=False) 

    if self.grow:
      olex.m('fuse')    

    atom.site=new_location_plus    
    olx.xf.au.SetAtomCrd(id, *new_location_plus)    # olx.xf.au.SetAtomU(id, *u_trans)
    olx.xf.EndUpdate()           

    if self.grow:
      olex.m('pack cell')
      olex.m('grow')    

    #this will hide the .hkl file in Wfn_job, so NoSpherA2 will not make a backup dir and tehrefore force ORCA to restart using the old wavefunction
    if keep_old_files == True:
      path = os.path.join(OV.FilePath(),"olex2","Wfn_job")
      try:
        shutil.move(os.path.join(path,olx.FileName()+".hkl"),os.path.join(path,olx.FileName()+".hidden_hkl"))
      except:
        pass  

    if hybrid_skip:
      f2_plus=f1_plus
      fcx2_plus=fcx1_plus
    else:      
      NoSpherA2.NoSpherA2.NoSpherA2_instance.launch() 
      f2_plus = self.get_fc_abs()   
      fcx2_plus=self.get_fc_complex() 


    if self.grow:
      olex.m('fuse')    

    atom.site=new_location_minus    
    olx.xf.au.SetAtomCrd(id, *new_location_minus)    
    olx.xf.EndUpdate() 

    if self.grow:
      olex.m('pack cell')
      olex.m('grow') 

    #this will hide the .hkl file in Wfn_job, so NoSpherA2 will not make a backup dir and tehrefore force ORCA to restart using the old wavefunction
    if keep_old_files == True:
      path = os.path.join(OV.FilePath(),"olex2","Wfn_job")
      try:
        shutil.move(os.path.join(path,olx.FileName()+".hkl"),os.path.join(path,olx.FileName()+".hidden_hkl"))
      except:
        pass      
    if hybrid_skip:
      f2_minus=f1_minus
      fcx2_minus=fcx1_minus
    else:
      NoSpherA2.NoSpherA2.NoSpherA2_instance.launch() 
      f2_minus = self.get_fc_abs()  
      fcx2_minus=self.get_fc_complex()     



    if self.grow:
      olex.m('fuse')        

    atom.site=old_location
    olx.xf.au.SetAtomCrd(id, *old_location)    
    olx.xf.EndUpdate()   


    if self.grow:
      olex.m('pack cell')
      olex.m('grow')     

    grad_cx_tru=(fcx2_plus-fcx2_minus)/(2*fullsize)
    grad_cx_approx=(fcx1_plus-fcx1_minus)/(2*fullsize) #with f0's tsc
    grad_cx_spher=(fcxsph_plus-fcxsph_minus)/(2*fullsize)   

    #2Re F* times df/fx

    grad_Y_tru=2*np.array((fcx0.conjugate()*grad_cx_tru).real)
    grad_Y_app=2*np.array((fcx0.conjugate()*grad_cx_approx).real)
    grad_Y_sph=2*np.array((fcxsph_0.conjugate()*grad_cx_spher).real)


    return (grad_Y_tru,grad_Y_app,grad_Y_sph)

    #this was the old version where the 'squared' gradients were calculated directly.

    grad_true = (f2_plus**2-f2_minus**2)/(2*fullsize) #with recalculated tsc
    grad_approx=(f1_plus**2-f1_minus**2)/(2*fullsize) #with f0's tsc
    grad_spher=(f_spher_plus**2-f_spher_minus**2)/(2*fullsize)

    return (grad_true,grad_approx, grad_spher)


  def get_fc_abs(self,table=True,original=False):
    #retrieves the |Fc| from olex2.
    fmr=FullMatrixRefine()
    #fmr=self.fmr
    if original:
      table_name= 'originaltsc.tsc'
      normal_eqns = fmr.run(build_only=True,table_file_name = table_name , normal_equations_class=leverage.leverage_normal_eqns)       
    elif table:
      molname=(os.path.split(OV.HKLSrc())[-1]).rsplit('.',1)[0]
      table_name=str(molname+".tsc")
      normal_eqns = fmr.run(build_only=True,table_file_name = table_name , normal_equations_class=leverage.leverage_normal_eqns)
    else:
      normal_eqns = fmr.run(build_only=True, normal_equations_class=leverage.leverage_normal_eqns)

    fcalcs=np.sqrt(np.array(fmr.get_fo_sq_fc(one_h_function=normal_eqns.one_h_linearisation)[1].as_intensity_array().data()))

    return fcalcs


  def get_fc_complex(self,table=True,original=False):
    #retrieves the complex fcalcs from olex2
    result=get_design(table,original)
    fcalcs2=[]

    return np.array(result.f_calc())

  def get_fo_fc_hkl(self,table=True,original=False):
    #retrieves Fo^2, Fc^2 and the hkl from olex2

    fmr = FullMatrixRefine()
    molname=(os.path.split(OV.HKLSrc())[-1]).rsplit('.',1)[0]
    table_name=str(molname+".tsc") #ZZZZ botch
    if original:
      table_name='originaltsc.tsc'
      normal_eqns = fmr.run(build_only=True,table_file_name = table_name , normal_equations_class=leverage.leverage_normal_eqns)
    elif table:
      normal_eqns = fmr.run(build_only=True,table_file_name = table_name , normal_equations_class=leverage.leverage_normal_eqns)
    else:
      normal_eqns = fmr.run(build_only=True, normal_equations_class=leverage.leverage_normal_eqns)

    intermediate=fmr.get_fo_sq_fc(one_h_function=normal_eqns.one_h_linearisation)
    fcalcs=np.array(intermediate[1].as_intensity_array().data())
    fobs=np.array(intermediate[0].as_intensity_array().data())
    hkl=np.array(intermediate[0].as_intensity_array().indices())

    return (fobs,fcalcs,hkl)  

  def refinement_step(self,hybrid=False):
    #performs one step of numerical refinement by calculating and applying a shift
    recordname,hklFile=self.recordname()

    At,Aa,Ah,As,Des_t,Des_s,Fc_t,W_t,Fc_s,W_s=retrieveData(recordname)    
    Fc_t=Fc_t[0]
    W_t=W_t[0]
    Fc_s=Fc_s[0]
    W_s=W_s[0]
    weightMatrix=np.diag(W_t)  
    fo_sq=self.get_fo_fc_hkl()[0]
    if hybrid:
      shift=getMyShift(Ah,fo_sq,Fc_t,W_t)
      uncertainties=getUncertainties(Ah,fo_sq,Fc_t,W_t)
    else:
      shift=getMyShift(At,fo_sq,Fc_t,W_t)
      uncertainties=getUncertainties(Ah,fo_sq,Fc_t,W_t)
    self.apply_shift(shift, uncertainties)


  def sph_ns_refine(self):
    #refines in the 'mixed' way using spherical derivatives and nonspherical Fcs. 
    #this performs almost as well as approximate nonspherical refinement, but takes the same time, so we do not see any practical applications
    shift,uncertainties=self.get_sph_ns_shift()
    self.apply_shift(shift,uncertainties)

  def apply_shift(self,shift,uncertainties):
    #applys a shift to a model in olex2
    #this is not an official way to do so but it works
    model = self.model  
    atoms = model.scatterers()    
    i=0

    if self.grow:
      olex.m('fuse')            

    for atomno in range(len(atoms)):
      dirshift=shift[i:i+3]
      atom=atoms[atomno]
      old_location=atom.site 
      new_location=old_location+dirshift  
      id = self.olx_atoms.atom_ids[atomno]  
      atom.site=new_location
      olx.xf.au.SetAtomCrd(id, *new_location)    
      olx.xf.EndUpdate()     
      i=i+3
      if atom.flags.use_u_iso():
        dirshift=shift[i]
        old_u= atom.u_iso 
        new_u=old_u+dirshift
        atom.u_iso=new_u
        olx.xf.au.SetAtomU(id, new_u)
        olx.xf.EndUpdate()          
        i=i+1
      if atom.flags.use_u_aniso():   
        dirshift=shift[i:i+6]
        old_u = atom.u_star
        ouc=adptbx.u_star_as_u_cart(self.xray_structure().unit_cell(), atoms[atomno].u_star)
        old_u_cart=(ouc[0],ouc[1],ouc[2],ouc[5],ouc[4],ouc[3])
        new_u=old_u+dirshift
        atom.u_star=new_u
        upc=adptbx.u_star_as_u_cart(self.xray_structure().unit_cell(), atoms[atomno].u_star)
        new_u_cart=(upc[0],upc[1],upc[2],upc[5],upc[4],upc[3])
        olx.xf.au.SetAtomU(id, *new_u_cart)
        olx.xf.EndUpdate()
        i=i+6 
    if self.refine_extinction and len(shift)-i==1: #this is a baaaad if. But will work for the current examples.
      exti = OV.GetExtinction()
      newexti=exti+shift[i]
      OV.SetExtinction(newexti,np.sqrt(uncertainties[i,i]))

    if self.grow:
      olex.m('pack cell')
      olex.m('grow')     


  def get_sph_ns_shift(self):
    #retrieves the shift for the mixed method which uses spherical derivatives and nonspherical Fcalcs
    #this performs almost as well as approximate nonspherical refinement, but takes the same time, so we do not see any practical applications    
    NoSpherA2.NoSpherA2.NoSpherA2_instance.launch()
    recordname,hklFile=self.recordname()
    fo_sq=self.get_fo_fc_hkl()[0]
    fc_sq=self.get_fc_abs()**2
    initial_design=get_design(table=False)
    index_count=(np.array(initial_design.f_calc())).size
    parameter_count=(np.array(initial_design.design_matrix())).size//index_count
    design=np.array(initial_design.design_matrix()).reshape(index_count,parameter_count)
    weights=np.array(initial_design.weights()) 
    shift=getMyShift(design, fo_sq, fc_sq, weights)
    uncertainties=getUncertainties(design,fo_sq,fc_sq,weights)   
    fc_sph=self.get_fc_abs(table=False)**2
    uncert_sph=getUncertainties(design, fo_sq, fc_sph, weights)
    return (shift,uncertainties)



  def app_shift(self,internal=False):
    #returns the approximate nonspherical shift vector and its uncertainties
    #there are two methods here, one using olex2's internal calculations and one using my own, to allow verification that it was the same.
    recordname,hklFile=self.recordname()
    fo_sq=self.get_fo_fc_hkl()[0]
    ##fc_sq=self.get_fc_abs()**2
    initial_design=get_design(table=True)
    weights=np.array(initial_design.weights()) 
    if internal:
      design=initial_design.design_matrix()
      fc_sq=self.get_fc_abs()**2
      shift=getOlexShift(design, fo_sq, fc_sq, weights)
      index_count=(np.array(initial_design.f_calc())).size
      parameter_count=(np.array(initial_design.design_matrix())).size//index_count
      design=np.array(initial_design.design_matrix()).reshape(index_count,parameter_count)
    else:
      design,fc=self.my_analytical_design_fcalc() 
      fc_sq=np.real(fc*np.conjugate(fc))
      shift=getMyShift(design, fo_sq, fc_sq, weights)  
    uncertainties=getUncertainties(design,fo_sq,fc_sq,weights)    
    return (shift,uncertainties)

  def my_analytical_design_fcalc(self):
    #a personal method to allow myself to calculate the design matrix, and thus to look at each step and troubleshoot. Essentially unused in the current format.
    hkls=self.get_fo_fc_hkl()[2]    
    parameters=self.retrieveParameters()
    formhkl,formfactors=retrieveFormFactors()
    atoms=self.retrieve_all_atoms()

    #this is a holdover from attempts to deal with extinction
    extinctionPresent=True
    fmr=FullMatrixRefine()
    extinction=fmr.run(build_only=True, normal_equations_class=leverage.leverage_normal_eqns).reparametrisation.extinction
    try:
      extinction=extinction.value
    except:
      extinction=0
      extinctionPresent=False


    hkllength=len(hkls)
    if extinctionPresent:
      paramlength=len(atoms)*9+1
    else:
      paramlength=len(atoms)*9


    theirfcalc=self.get_fc_complex()

    fo2 = self.reflections.f_sq_obs_filtered
    symmetry_ops=fo2.crystal_symmetry().space_group().all_ops()        
    symmetry_matrices=[]
    symmetry_shifts=[]
    for matrix in symmetry_ops:
      symmetry_matrices+=[np.reshape(matrix.as_double_array()[0:9], (3,3))]
      symmetry_shifts+=[np.reshape(matrix.as_double_array()[9:], (3))]

    derivative=0
    fcalcs=np.full((hkllength),np.nan,dtype='complex128')
    adjusted_fcalcs=np.full((hkllength),np.nan,dtype='complex128')
    design=np.full((hkllength,paramlength),np.nan,dtype='complex128')
    design_exti_adjusted=np.full((hkllength,paramlength),np.nan,dtype='complex128')
    design_exti_from_intensity=np.full((hkllength,paramlength),np.nan,dtype='complex128')
    design_with_their_fcalc=np.full((hkllength,paramlength),np.nan,dtype='complex128')
    design_with_their_fcalc_noim=np.full((hkllength,paramlength),np.nan,dtype='complex128')
    design_with_their_fcalc_nonconj=np.full((hkllength,paramlength),np.nan,dtype='complex128')
    design_with_their_fcalc_fully=np.full((hkllength,paramlength),np.nan,dtype='complex128')
    sin2theta=np.full(hkllength,np.nan)


    orthogonalization=np.reshape(np.array(self.model.unit_cell().orthogonalization_matrix()),(3,3))
    orthogonalization_inv=np.linalg.inv(orthogonalization) 
    wavelength=self.wavelength


    def sintwotheta(h):
      dstar=np.dot(orthogonalization_inv.T,h)
      dstarsquare=np.dot(dstar,dstar)
      x=dstarsquare*wavelength**2/4 #sinsquare
      sintwotheta=2*np.sqrt(x*(1-x))
      return sintwotheta

    extifcmultipliers=np.full(hkllength, np.nan)
    extiderivmultipliers=np.full(hkllength, np.nan)

    for k, hkl in enumerate(hkls):
      fcalcmycalc=complex(0,0)
      sin2theta[k]=sintwotheta(hkl)
      for i, atom in enumerate(atoms):
        formfac=formfactors[k][i]+complex(atom.fp,atom.fdp)
        derivative=np.full((9),complex(0.,0.))
        choicexyz=parameters[i*9:i*9+3]
        choiceu=parameters[i*9+3:i*9+9]
        for j, matrix in enumerate(symmetry_matrices):
          rotatedh=np.dot(matrix,hkl).astype(int)
          try:
            rotatedhloc=np.intersect1d(np.where(formhkl[:,0]==rotatedh.astype(int)[0]),np.intersect1d(np.where(formhkl[:,1]==rotatedh.astype(int)[1]),np.where(formhkl[:,2]==rotatedh.astype(int)[2])))[0]     
            formfac=formfactors[rotatedhloc][i]+complex(atom.fp,atom.fdp)
          except:
            rotatedh=-rotatedh
            rotatedhloc=np.intersect1d(np.where(formhkl[:,0]==rotatedh.astype(int)[0]),np.intersect1d(np.where(formhkl[:,1]==rotatedh.astype(int)[1]),np.where(formhkl[:,2]==rotatedh.astype(int)[2])))[0]     
            formfac=np.conjugate(formfactors[rotatedhloc][i])+complex(atom.fp,atom.fdp)
            rotatedh=-rotatedh    
          fcalcContrib=formfac*np.exp(-2*np.pi**2*np.dot(np.array(rotatedh),np.dot(np.array([[choiceu[0],choiceu[3],choiceu[4]],[choiceu[3],choiceu[1],choiceu[5] ],[choiceu[4],choiceu[5],choiceu[2] ]]),np.array(rotatedh))))*np.exp(2*np.pi*complex(0,1)*np.dot(rotatedh,choicexyz))*np.exp(complex(0,1)*2*np.pi*np.dot(rotatedh,symmetry_shifts[j]))
          derivative[0]+=fcalcContrib*complex(0,1)*2*np.pi*rotatedh[0]
          derivative[1]+=fcalcContrib*complex(0,1)*2*np.pi*rotatedh[1]
          derivative[2]+=fcalcContrib*complex(0,1)*2*np.pi*rotatedh[2]
          derivative[3]+=-2*np.pi**2*fcalcContrib*rotatedh[0]**2
          derivative[4]+=-2*np.pi**2*fcalcContrib*rotatedh[1]**2
          derivative[5]+=-2*np.pi**2*fcalcContrib*rotatedh[2]**2
          derivative[6]+=-4*np.pi**2*fcalcContrib*rotatedh[0]*rotatedh[1]
          derivative[7]+=-4*np.pi**2*fcalcContrib*rotatedh[0]*rotatedh[2]
          derivative[8]+=-4*np.pi**2*fcalcContrib*rotatedh[1]*rotatedh[2]          
          choicexyz=parameters[i*9:i*9+3]
          choiceu=parameters[i*9+3:i*9+9]
          fcalcmycalc+=fcalcContrib
        design[k,i*9:(i+1)*9]=derivative
      design_with_their_fcalc[k]=design[k]*2*np.conjugate(theirfcalc[k])
      exti_mult_factor=1+0.001*wavelength**3*extinction*(fcalcmycalc*np.conjugate(fcalcmycalc))/sin2theta[k]
      extifcmultipliers[k]=np.power(exti_mult_factor,-0.5)
      extiderivmultipliers[k]=np.power(exti_mult_factor,-0.5)/2*(1+1/exti_mult_factor)
      design_exti_adjusted[k]=np.power(exti_mult_factor,-0.25)*(design[k]-(fcalcmycalc*0.001*wavelength**3*extinction*2*np.real(fcalcmycalc.conjugate()*design[k]))/(sin2theta[k]*4*exti_mult_factor))
      design_with_their_fcalc_nonconj[k]=design[k]*2*theirfcalc[k]#np.conjugate(theirfcalc[k])
      design_with_their_fcalc_noim=np.real(design[k])*2*np.real(theirfcalc[k])
      design[k]*=2*np.conjugate(fcalcmycalc)
      design[k]=np.real(design[k])
      fcalcs[k]=fcalcmycalc
      adjusted_fcalcs[k]=fcalcmycalc*np.power(exti_mult_factor,-0.25 )
      design_exti_adjusted[k]*=2*np.conjugate(adjusted_fcalcs[k])
      design_exti_from_intensity[k]=design[k]*np.power(exti_mult_factor,-0.5)*(1-(fcalcmycalc*np.conjugate(fcalcmycalc)*0.001*extinction*wavelength**3)/(sin2theta[k]*2*exti_mult_factor))
      if extinctionPresent:
        design[k,paramlength-1]=-np.power(exti_mult_factor,-1.5)*(fcalcmycalc*np.conjugate(fcalcmycalc))**2*0.001*wavelength**3/(2*sin2theta[k])
        design_exti_from_intensity[k,paramlength-1]=design[k,paramlength-1]
        design_exti_adjusted[k,paramlength-1]=design[k,paramlength-1]
    design=np.real(design)
    design_exti_adjusted=np.real(design_exti_adjusted)
    design_with_their_fcalc=np.real(design_with_their_fcalc) ##ZZZZ 
    design_exti_from_intensity=np.real(design_exti_from_intensity)
    adjusted_fcalc=fcalcs*np.power(1+0.001*wavelength**3*extinction*(fcalcs*np.conjugate(fcalcs))/sin2theta,-0.25 )
    return(design,adjusted_fcalc)



  def next_shift(self,hybrid=False):
    #finds and returns the next numerical shift  vector

    molname=(os.path.split(OV.HKLSrc())[-1]).rsplit('.',1)[0]
    recordname=molname+"_gradient_records.txt"


    At,Aa,Ah,As,Des_t,Des_s,Fc_t,W_t,Fc_s,W_s=retrieveData(recordname)    
    Fc_t=Fc_t[0]
    W_t=W_t[0]
    Fc_s=Fc_s[0]
    W_s=W_s[0]
    weightMatrix=np.diag(W_t)  

    fo_sq=self.get_fo_fc_hkl()[0]

    if hybrid:
      shift=getMyShift(Ah,fo_sq,Fc_t,W_t)
    else:
      shift=getMyShift(At,fo_sq,Fc_t,W_t)
    return shift
  def getAllShifts(self):
    #finds and returns the next shift vectors for numerical, hybrid, approximate and spherical refinement, their uncertainties, and the shifts and uncertainties for mixed reifnement
    molname=(os.path.split(OV.HKLSrc())[-1]).rsplit('.',1)[0]
    recordname=molname+"_gradient_records.txt"


    At,Aa,Ah,As,Des_t,Des_s,Fc_t,W_t,Fc_s,W_s=retrieveData(recordname)    
    Fc_t=Fc_t[0]
    W_t=W_t[0]
    Fc_s=Fc_s[0]
    W_s=W_s[0]
    weightMatrix=np.diag(W_t)  

    fo_sq=self.get_fo_fc_hkl()[0]


    shiftH=getMyShift(Ah,fo_sq,Fc_t,W_t)
    uncertaintiesH=getUncertainties(Ah,fo_sq,Fc_t,W_t,expanded=True) 
    shiftT=getMyShift(At,fo_sq,Fc_t,W_t)
    uncertaintiesT=getUncertainties(At,fo_sq,Fc_t,W_t,expanded=True) 
    shiftA=getMyShift(Aa,fo_sq,Fc_t,W_t)
    uncertaintiesA=getUncertainties(Aa,fo_sq,Fc_t,W_t,expanded=True) 
    shiftS=getMyShift(As,fo_sq,Fc_s,W_s)
    uncertaintiesS=getUncertainties(As,fo_sq,Fc_s,W_s,expanded=True) 
    shiftM=getMyShift(As,fo_sq,Fc_t,W_t)
    uncertaintiesM=getUncertainties(As,fo_sq,Fc_t,W_t,expanded=True) 

    return (shiftT,shiftH,shiftA,shiftS,uncertaintiesT,uncertaintiesH,uncertaintiesA,uncertaintiesS,shiftM,uncertaintiesM)    
  def get_full_data(self,steps=0, hybrid=False,mix=False):
    #calls the other function ('get_data_each_shift' below), and puts in atomic data into the comparisonData file
    self.get_data_each_shift(step=steps,hybrid=hybrid,mix=mix)
    model = self.model  
    atoms = model.scatterers()  

    for atom in atoms:
      atomParams=[atom.label,atom.site]
      if atom.flags.use_u_iso():
        atomParams+=[atom.u_iso]
      elif atom.flags.use_u_aniso():
        atomParams+=[atom.u_star]
      depositFile=open(self.depositFile,'a')
      depositFile.write(str(atomParams)) #ZZZZ EXTI write out here?
      depositFile.write("\n")
      depositFile.close()

  def retrieveParameters(self):
    #gets a list of parameters for all atoms
    model = self.model  
    atoms = model.scatterers()  
    paramlist=[]
    for atom in atoms:
      atomParams=list(atom.site)
      if atom.flags.use_u_iso():
        atomParams+=[atom.u_iso]
      elif atom.flags.use_u_aniso():
        atomParams+=atom.u_star
      paramlist+=atomParams

    return paramlist    
  def retrieve_all_atoms(self):
    #returns a list of all atoms in the model.
    model = self.model  
    atoms = model.scatterers()  
    return atoms

  def get_data_each_shift(self, step=0,hybrid=False,mix=False):
    #retrieves and stores a lot of information in comparisonData such as fcalcs, maximal shifts/esd
    model = self.model  
    atoms = model.scatterers()
    shifts=self.getAllShifts()    
    metrical=self.metrical

    recordname,hklFile=self.recordname()
    At,Aa,Ah,As,Des_t,Des_s,Fc_t,W_t,Fc_s,W_s=retrieveData(recordname)    
    Fc_t=Fc_t[0]
    W_t=W_t[0]
    Fc_s=Fc_s[0]
    W_s=W_s[0] 
    fo_sq=self.get_fo_fc_hkl()[0]
    sigma=np.sqrt(1/np.array(W_t))
    uncertainties=getUncertainties(At, fo_sq, Fc_t, W_t)
    depositFile=open(self.depositFile,'a')
    shiftlist=self.getAllShifts()
    shift=shiftlist[0]
    uncert=shiftlist[4][1]
    if hybrid:
      shift=shiftlist[1]
      uncert=shiftlist[5][1]
    if mix:
      shift=shiftlist[8]
      uncert=shiftlist[9][1]
    shiftoveruncert=shift/uncert    


    #We are entering Big Botch Zone.
    height=shift.shape[0]//9
    reshapeshift=np.reshape(shift, (height,9))
    reshapese=np.reshape(shiftoveruncert,(height,9))
    #maximal values for shifts and shift/esd
    maxse=np.max(abs(shiftoveruncert))
    maxpos=np.max(abs(reshapeshift[:,:3]))
    maxadp=np.max(abs(reshapeshift[:,3:]))
    locse=np.where(abs(reshapese)==maxse)
    locpos=np.where(abs(reshapeshift)==maxpos)
    locadp=np.where(abs(reshapeshift)==maxadp)
    labels=["x","y","z","u11","u22","u33","u12","u13","u23"]
    labelse=labels[locse[1][0]]
    labelpos=labels[locpos[1][0]]
    labeladp=labels[locadp[1][0]]


    depositFile.write("Step "+str(step)+": Max Shift="+str(maxpos)+" [Pos] ("+atoms[int(locpos[0][0])].label+"."+labelpos+"), "+str(maxadp)+" [ADP] ("+atoms[int(locadp[0][0])].label+"."+ labeladp+"), Shift/esd:"+ str(maxse)+" ("+atoms[int(locse[0][0])].label+"."+labelse+")")
    depositFile.write("\n")
    depositFile.write(str(["R factors: ", rfac(Fc_t, fo_sq, W_t, sigma)]))
    depositFile.write(str(["R factors without cutoff: ", rfac(Fc_t, fo_sq, W_t, sigma,thresh=0)]))
    depositFile.write(str(["R factor Spher: ", rfac(Fc_s, fo_sq, W_s, sigma)]))
    depositFile.write("\n")

    depositFile.close()

  def get_rfac_data(self):
    #generates the rfactors for the current structure and stores them in comparisonData
    recordname,hklFile=self.recordname()
    At,Aa,Ah,As,Des_t,Des_s,Fc_t,W_t,Fc_s,W_s=retrieveData(recordname)    
    Fc_t=Fc_t[0]
    W_t=W_t[0]
    Fc_s=Fc_s[0] 
    W_s=W_s[0]
    fo_sq=self.get_fo_fc_hkl()[0]
    sigma=np.sqrt(1/np.array(W_t))


    fo_sq=self.get_fo_fc_hkl()[0]
    Fc_t=self.get_fc_abs()**2
    Fc_s=self.get_fc_abs(table=False)
    initial_design=get_design(table=False)
    index_count=(np.array(initial_design.f_calc())).size
    parameter_count=(np.array(initial_design.design_matrix())).size//index_count
    design=np.array(initial_design.design_matrix()).reshape(index_count,parameter_count)
    W_t=np.array(initial_design.weights()) 
    sigma=np.sqrt(1/np.array(W_t))    


    depositFile=open(self.depositFile,'a')
    depositFile.write(str(["R factors: ", rfac(Fc_t, fo_sq, W_t, sigma)]))
    depositFile.write(str(["R factors without cutoff: ", rfac(Fc_t, fo_sq, W_t, sigma,thresh=0)]))
    depositFile.write(str(["R factor Spher: ", rfac(Fc_s, fo_sq, W_s, sigma)]))
    depositFile.write("\n")    


  def get_rfac_data_mix(self):
    #gets the rfactors (into comparisonData) without requiring the generation of a gradient file.
    NoSpherA2.NoSpherA2.NoSpherA2_instance.launch()
    recordname,hklFile=self.recordname()
    fo_sq=self.get_fo_fc_hkl()[0]
    Fc_t=self.get_fc_abs()**2
    Fc_s=self.get_fc_abs(table=False)
    initial_design=get_design(table=False)
    index_count=(np.array(initial_design.f_calc())).size
    parameter_count=(np.array(initial_design.design_matrix())).size//index_count
    design=np.array(initial_design.design_matrix()).reshape(index_count,parameter_count)
    W_t=np.array(initial_design.weights()) 
    sigma=np.sqrt(1/np.array(W_t))    

    depositFile=open(self.depositFile,'a')
    depositFile.write(str(["R factors: ", rfac(Fc_t, fo_sq, W_t, sigma)]))
    depositFile.write(str(["R factors without cutoff: ", rfac(Fc_t, fo_sq, W_t, sigma,thresh=0)]))
    depositFile.write(str(["R factor Spher: ", rfac(Fc_s, fo_sq, W_t, sigma)]))
    depositFile.write("\n")        


  def get_base_data(self):
    #records initial data/settings into comparisonData, such as the options used in NoSpherA2 and what epsilon is being used.
    depositFile=open(self.depositFile,'a')
    depositFile.write("New Record\n")
    depositFile.write("Orthogonalization Matrix\n")
    depositFile.write(str(self.orthogonalization))
    depositFile.write("\n")    
    #OV.GetParam('snum.NoSpherA2.use_aspherical')#,True)
    depositFile.write("Source: "+str(OV.GetParam('snum.NoSpherA2.source'))+"\n") #,"ORCA")
    depositFile.write("Calculate: "+str(OV.GetParam('snum.NoSpherA2.Calculate'))+"\n")#,True)
    depositFile.write("Basis Name: "+str(OV.GetParam('snum.NoSpherA2.basis_name'))+"\n")#,"def2-TZVPP")
    depositFile.write("Becke Accuracy: "+str(OV.GetParam('snum.NoSpherA2.becke_accuracy'))+"\n")#,"Normal")
    depositFile.write("Method: "+str(OV.GetParam('snum.NoSpherA2.method'))+"\n")#,"PBE")
    depositFile.write("Full HAR: "+str(OV.GetParam('snum.NoSpherA2.full_HAR'))+"\n")#,True)
    depositFile.write("Epsilon: "+str(OV.GetParam('lauranonsphertests.vars.epsilon'))+"\n")#,True)
    depositFile.close()

  def deposit_phrase(self, phrase):
    #puts a phrase into comparisonData
    depositFile=open(self.depositFile,'a')
    depositFile.write(phrase)
    depositFile.write("\n")
    depositFile.close()


  def recordname(self):
    #returns the record name associated to the gradient records for the structure, and the hkl file I use 
    #(as sometimes the hkl file is out of order with respect to the order olex2 gives the data in, so one must make a rearranged duplicate)
    molname=(os.path.split(OV.HKLSrc())[-1]).rsplit('.',1)[0]
    recordname=molname+"_gradient_records.txt"
    hklFile=molname+"2.hkl"
    try:
      x=open(hklFile,'r')
    except:
      hklFile=molname+".hkl"  
    return (recordname,hklFile)
  def compare_structures(self,alternate):
    #compares two structures with regards to their wR2 and returns the structure with the lowest wR2. This is for the purpose of returning to the 'best' structure.
    recordname,hklFile=self.recordname()
    At,Aa,Ah,As,Des_t,Des_s,Fc_t,W_t,Fc_s,W_s=retrieveData(recordname)    
    Fc_t=Fc_t[0]
    W_t=W_t[0]
    Fc_s=Fc_s[0]
    W_s=W_s[0]
    fo_sq=self.get_fo_fc_hkl()[0]
    sigma=np.sqrt(1/np.array(W_t))

    fo_sq=self.get_fo_fc_hkl()[0]

    newwR2=rfac(Fc_t, fo_sq, W_t, sigma,thresh=0)[1] 
    wR2alt=alternate.wR2
    newatoms=copy.deepcopy(self.model.scatterers())

    if newwR2<wR2alt:
      replacement = storedStructure(newatoms,newwR2)
      return (replacement,True)
    return (alternate,False)
  def replace_structure(self,newStructure):
    #called at the end of a refinement loop, replaces the parameters of all atoms with those associated to newStructure - usually the one with lowest wR2.
    newAtoms=newStructure.atoms


    model = self.model  
    atoms = model.scatterers()    
    i=0 
    for atomno in range(len(atoms)):
      atom=atoms[atomno]
      newatom=newAtoms[atomno]
      old_location=atom.site 
      new_location=newatom.site
      id = self.olx_atoms.atom_ids[atomno]  
      atom.site=new_location
      olx.xf.au.SetAtomCrd(id, *new_location)    
      olx.xf.EndUpdate()     
      i=i+3
      if atom.flags.use_u_iso():
        old_u= atom.u_iso 
        new_u=newatom.u_iso
        atom.u_iso=new_u
        olx.xf.au.SetAtomU(id, new_u)
        olx.xf.EndUpdate()          
        i=i+1
      if atom.flags.use_u_aniso():   
        old_u = atom.u_star
        new_u=newatom.u_star


        ouc=adptbx.u_star_as_u_cart(self.xray_structure().unit_cell(), atoms[atomno].u_star)
        old_u_cart=(ouc[0],ouc[1],ouc[2],ouc[5],ouc[4],ouc[3])
        atom.u_star=new_u
        upc=adptbx.u_star_as_u_cart(self.xray_structure().unit_cell(), atoms[atomno].u_star)
        new_u_cart=(upc[0],upc[1],upc[2],upc[5],upc[4],upc[3])
        olx.xf.au.SetAtomU(id, *new_u_cart)
        olx.xf.EndUpdate()
        i=i+6  
  def make_own_hkl(self,fcalcs=False):
    #prints a hkl-ish formatted list in the same order as stored in olex2, which may not match that in the associated .hkl file.
    fo_sq,fc_sq,hkls=self.get_fo_fc_hkl(table=False)
    #fo_sq=self.get_fo_fc_hkl(table=False)[0]
    #hkls=self.get_fo_fc_hkl(table=False)[2]
    initial_design=get_design(table=False)    
    weights=initial_design.weights()
    sigma=np.sqrt(1/np.array(weights))
    if fcalcs:
      fo_sq=fc_sq
      sigma=fc_sq*0.05
    for i,line in enumerate(hkls):
      foformat='{:> .6g}'
      sigformat='{:> .6g}'
      if fo_sq[i]<1:
        foformat='{:> .5f}'
      if sigma[i]<1:
        sigformat='{:> .5f}'
  
      print( '{:4d}'.format(int(line[0]))+'{:4d}'.format(int(line[1]))+'{:4d}'.format(int(line[2]))+foformat.format(fo_sq[i]).rjust(8)+sigformat.format(sigma[i]).rjust(8))
    print("   0   0   0 ")
    
      
    #for i,line in enumerate(hkls):
      #print(str(line[0])+" "+str(line[1])+" "+str(line[2])+" "+str(fo_sq[i])+" "+str(sigma[i]))



def get_design(table=True,original=False):
  #pulls the design matrix from olex2
  #this is very much directly copied from internal files, I cannot explain how it works!
  fmr = FullMatrixRefine()
  molname=(os.path.split(OV.HKLSrc())[-1]).rsplit('.',1)[0]
  table_name=str(molname+".tsc") #ZZZZ botch
  if original:
    table_name='originaltsc.tsc'
    normal_eqns = fmr.run(build_only=True,table_file_name = table_name , normal_equations_class=leverage.leverage_normal_eqns)
  elif table:
    normal_eqns = fmr.run(build_only=True,table_file_name = table_name , normal_equations_class=leverage.leverage_normal_eqns)
  else:
    normal_eqns = fmr.run(build_only=True, normal_equations_class=leverage.leverage_normal_eqns)

  if normal_eqns.f_mask is not None:
    f_mask = normal_eqns.f_mask.data()
  else:
    f_mask = flex.complex_double()   

  extinction_correction = normal_eqns.reparametrisation.extinction
  if extinction_correction is None:
    extinction_correction = xray.dummy_extinction_correction()    

  def args(scale_factor, weighting_scheme):
    args = (normal_eqns,
            normal_eqns.observations,
          f_mask,
          weighting_scheme,
          scale_factor,
          normal_eqns.one_h_linearisation,
          normal_eqns.reparametrisation.jacobian_transpose_matching_grad_fc(),
          extinction_correction)
    return args

  normal_eqns.reparametrisation.linearise()
  normal_eqns.reparametrisation.store()
  scale_factor = float(olx.xf.rm.OSF())
  scale_factor *= scale_factor
  result = ext.build_design_matrix(*args(scale_factor,normal_eqns.weighting_scheme))  

  return result




def parameter_labels(fmr, n_params):
  #retrieves the labels for each parameter in the structure vector x
  #I don't think we ever actually deal with the structure vector, but it points out which parameters are independant of one another.
  #Jt is the jacobean-transpose matrix which takes x to y.
  annotations = [x for x in fmr.reparametrisation.component_annotations]
  annotations_1 = []
  labels = []
  ann_1_idx = 0
  if fmr.reparametrisation.twin_fractions is not None:
    basf_n = 1
    for fraction in fmr.reparametrisation.twin_fractions:
      if fraction.grad:
        annotations_1.append("BASF%s" %basf_n)
        basf_n += 1
  if fmr.reparametrisation.extinction is not None and fmr.reparametrisation.extinction.grad:
    annotations_1.append("EXTI")
  Jt = fmr.reparametrisation.jacobian_transpose_matching_grad_fc()
  for j in range(n_params):
    label = []
    for k in range(Jt.n_cols):
      if Jt[(j,k)]:
        label.append("%s" %(annotations[k]))
    if len(label) == 0:
      label.append(annotations_1[ann_1_idx])
      ann_1_idx += 1
    labels.append(", ".join(label))
  return labels


def retrieveData(filename):
  #retrieves data from a gradient_records file
  test=open(filename,'r')

  grads_aspher_tru=[]
  grads_aspher_app=[]
  tsc_design_ma=[]
  grads_aspher_hybrid=[]
  grads_spher=[]
  design_ma=[]
  fcalcs_tru=[]
  weights_tru=[]
  fcalcs_spher=[]
  weights_spher=[]

  lines = test.readlines()

  i=0
  for line in lines:
    if line[0]!='[':
      i+=1
      continue

    arrayline=np.fromstring((line.rsplit("]")[0]).rsplit("[")[1],sep=',')

    if i==2: 
      grads_aspher_tru+= [list(arrayline)]
    elif i==4:
      grads_aspher_app+= [list(arrayline)]
    elif i==6:
      tsc_design_ma+= [list(arrayline)]
    elif i==8:
      grads_aspher_hybrid+= [list(arrayline)]
    elif i==10:
      grads_spher+= [list(arrayline)]
    elif i==12:
      design_ma+= [list(arrayline)]
    elif i==14:
      fcalcs_tru=[list(arrayline)]
    elif i==16:
      weights_tru=[list(arrayline)]
    elif i==18:
      fcalcs_spher=[list(arrayline)]
    elif i==20:
      weights_spher=[list(arrayline)]



  test.close()        

  grads_aspher_tru=np.array(grads_aspher_tru)
  grads_aspher_app=np.array(grads_aspher_app)
  grads_aspher_hybrid=np.array(grads_aspher_hybrid)
  grads_spher=np.array(grads_spher)
  tsc_design_ma=np.array(tsc_design_ma)
  design_ma=np.array(design_ma) 




  return (grads_aspher_tru,grads_aspher_app,grads_aspher_hybrid,grads_spher,tsc_design_ma,design_ma,fcalcs_tru,weights_tru,fcalcs_spher,weights_spher)    

def retrieveFormFactors():
  #pulls form factors from the .tsc file and returns them along with their hkl for verification.
  #This will not retrieve the atom names, and simply expects them to match the order olex2 gives them in.
  molname=(os.path.split(OV.HKLSrc())[-1]).rsplit('.',1)[0]
  filename=str(molname+".tsc")
  test=open(filename,'r')
  hkl=[]
  formfacs=[]
  for line in test:
    if line[0] in ['T','S','D']:
      continue
    values=np.array(line.split())
    hkl+=[values[:3]]
    formfacline=[]
    for entry in values[3:]:
      toTwo=np.fromstring(entry,sep=",")
      formfac=complex(toTwo[0],toTwo[1])
      formfacline+=[formfac]
    formfacs+=[formfacline]
  hkl=np.array(hkl).astype(int)
  return (hkl, formfacs)


def getMyShift(design, fo_sq, fc_sq, weights):
  #directly creates a shift vector from a design matrix, Fo^2, Fc^2 and a weighting scheme
  fmr=FullMatrixRefine()
  initial_fmr=fmr.run(build_only=True)
  reparam=fmr.reparametrisation
  reparam.linearise()
  jt=reparam.jacobian_transpose_matching_grad_fc()
  jt=jt.as_dense_matrix().as_numpy_array()  

  design=np.array(design)
  fo_sq=np.array(fo_sq)
  fc_sq=np.array(fc_sq)
  weights=np.array(weights)
  fcwfc=np.dot(weights*fc_sq,fc_sq)
  ktilde=np.dot(weights*fo_sq,fc_sq)/fcwfc
  dktilde=np.dot(design.T,weights*(fo_sq-2*ktilde*fc_sq))/fcwfc
  adjusted_design=ktilde*design+np.outer(fc_sq,dktilde)
  weights=np.diag(weights)       
  fdiff=fo_sq-ktilde*fc_sq     
  adw=np.dot(adjusted_design.T,weights)
  varianceMatrix=np.linalg.inv(np.dot(adw,adjusted_design))
  uncertainties=np.sqrt(np.diag(varianceMatrix)) #ZZZ uncertainties
  matrix=np.dot(np.linalg.inv(np.dot(adw,adjusted_design)),adw)

  shifts=np.dot(matrix,fdiff)   
  shifts_y=np.dot(jt.T,shifts)

  return(shifts_y)  

def getOlexShift(design, fo_sq, fc_sq, weights):
  #creates a shift vector using olex2's internal functions.
  from scitbx.lstbx import normal_eqns as normal_eqns_m
  fmr=FullMatrixRefine()
  initial_fmr=fmr.run(build_only=True)
  reparam=fmr.reparametrisation
  reparam.linearise()
  jt=reparam.jacobian_transpose_matching_grad_fc()
  jt=jt.as_dense_matrix().as_numpy_array()    

  fo_sq=flex.double(list(fo_sq))
  fc_sq=flex.double(list(fc_sq))
  weights=flex.double(list(weights))
  fmr = FullMatrixRefine()
  fmrrun = fmr.run(build_only=True, normal_equations_class=leverage.leverage_normal_eqns)
  test_ne = normal_eqns_m.non_linear_ls_with_separable_scale_factor_BLAS_3(fmrrun.reparametrisation.n_independents) #63 for epoxide
  test_ne.add_equations(fc_sq,design, fo_sq, weights)
  test_ne.finalise()
  test_ne.solve()
  vector=[x for x in test_ne.step()]
  shifts_y=np.dot(jt.T,vector)
  return vector    

def getUncertainties(design, fo_sq, fc_sq, weights,expanded=False,sigmas=None,thresh=2):
  #creates the variance-covariance matrix associated to the current minima/next shift
  fmr=FullMatrixRefine()
  initial_fmr=fmr.run(build_only=True)
  reparam=fmr.reparametrisation
  reparam.linearise()
  jt=reparam.jacobian_transpose_matching_grad_fc()
  jt=jt.as_dense_matrix().as_numpy_array()    

  design=np.array(design)
  fo_sq=np.array(fo_sq)
  fc_sq=np.array(fc_sq)
  weights=np.array(weights)


  fo_sq2=fo_sq
  fc_sq2=fc_sq
  weights2=weights

  fcwfc=np.dot(weights*fc_sq,fc_sq)
  ktilde=np.dot(weights*fo_sq,fc_sq)/fcwfc
  dktilde=np.dot(design.T,weights*(fo_sq-2*ktilde*fc_sq))/fcwfc
  adjusted_design=ktilde*design+np.outer(fc_sq,dktilde)
  multiplier = np.dot(weights*(ktilde*fc_sq-fo_sq),(ktilde*fc_sq-fo_sq))/(int(max(design.shape))-int(min(design.shape)))
  weights=np.diag(weights)      
  adw=np.dot(adjusted_design.T,weights)
  varianceMatrix=np.linalg.inv(np.dot(adw,adjusted_design))
  uncertainties=np.sqrt(np.diag(varianceMatrix)) #ZZZ uncertainties

  uncertainties=np.sqrt(multiplier)*uncertainties #ZZZZ #http://pd.chem.ucl.ac.uk/pdnn/refine1/errors.htm & http://pd.chem.ucl.ac.uk/pdnn/refine1/rfacs.htm
  adjVarMat=multiplier*varianceMatrix 
  adjVarMat_y=np.dot(jt.T,adjVarMat)
  uncert_y=np.dot(jt.T,uncertainties)
  if expanded:
    return (adjVarMat_y,uncert_y)
  return(adjVarMat_y)  


def rfac(fc_sq,fo_sq,w,sigmas,thresh=2):
  #gets the r1 and wr2 factor for the given data
  fo_sq=np.array(fo_sq)
  fc_sq=np.array(fc_sq)
  w=np.array(w)
  cutoff = thresh*np.array(sigmas)

  selection=np.where(fo_sq>=cutoff)
  fo_sq=fo_sq[selection]
  fc_sq=fc_sq[selection]
  w=w[selection]

  fcwfc=np.dot(w*fc_sq,fc_sq)
  ktilde=np.dot(w*fo_sq,fc_sq)/fcwfc   
  r1=np.sum(abs(abs(np.sqrt(fo_sq))-abs(np.sqrt(ktilde*fc_sq))))/np.sum(abs(np.sqrt(fo_sq)))
  r3=np.sum(abs(np.dot(w*(fo_sq-ktilde*fc_sq),fo_sq-ktilde*fc_sq)))/np.sum(abs(w*fo_sq**2))
  wr2=np.sqrt(r3) ##this is the wR2

  return r1, wr2

def run_many_tests():
  #runs nonspherical tests for different bases,  becke accuracies and methods
  #this hasn't actually been run
  import os
  OV.SetParam('snum.NoSpherA2.use_aspherical',True)
  OV.SetParam('snum.NoSpherA2.full_HAR',False)  
  #OV.SetParam('snum.NoSpherA2.source',"ORCA") #ZZZZ Removed set ORCA
  OV.SetParam('snum.NoSpherA2.Calculate',True)
  olex.m('anis -h')


  bases=["def2-TZVPP","def2-TZVPP"]
  methods=["PBE","HF","B3LYP"]
  beckuracy=["Normal","High","Low"] #it's a portmanteau of becke and accuracy hahaha I am funny

  for acc in beckuracy:
    OV.SetParam('snum.NoSpherA2.becke_accuracy',acc)
    for method in methods:
      OV.SetParam('snum.NoSpherA2.method',method)
      for basis in bases:
        OV.SetParam('snum.NoSpherA2.basis_name',basis)
        for j in range(20):
          olex.m('refine')    #refine to aspherical optimum  
          OV.SetParam('snum.NoSpherA2.Calculate',True)               
        olex.m('spy.LauraNonSpherTests.refine_multiple()')
        os.rename('comparisonData.txt','comparisonData_'+acc+'_'+method+'_'+basis+'.txt')
  shutil.rmtree(r"olex2\Wfn_job")
  shutil.rmtree(r"olex2\NoSpherA2_history")

def run_test_parameter_choices():
  #runs full numerical refinements for 20 steps for each epsilon in the list below (1e-1 to 1e-6 but in order of 'likely usefulness')
  import time
  tic = time.time()
  #first: disable NoSpherA2 and refine a few times to make sure things are settled
  OV.SetParam('snum.NoSpherA2.use_aspherical',False)
  olex.m('anis')
  olex.m('anis -h')
  for i in range(100):
    olex.m('refine')
  
  epsilonList=[1e-3,1e-4,1e-2,1e-5,1e-1,0.2,1e-6]
  beckuracyList=["Low","Normal"]#["High","Low","Normal"]
  basesList=["def2-TZVP","def2-SVP","cc-pVTZ"]
  

  #put NoSpherA2 settings
  OV.SetParam('snum.NoSpherA2.use_aspherical',True)
  OV.SetParam('snum.NoSpherA2.full_HAR',False)  
  OV.SetParam('snum.NoSpherA2.source',"DISCAMB")
  OV.SetParam('snum.NoSpherA2.Calculate',True)  
  #OV.SetParam('snum.NoSpherA2.basis_name',"def2-TZVPP")
  OV.SetParam('snum.NoSpherA2.becke_accuracy',"High")
  OV.SetParam('snum.NoSpherA2.method',"PBE")

  olex.m('spy.LauraNonSpherTests.make_gradient_file()')
  olex.m('spy.LauraNonSpherTests.deposit_phrase("Spherical Model")')
  olex.m('spy.LauraNonSpherTests.do_test()') 
   

  for beckuracy in beckuracyList:
    OV.SetParam('snum.NoSpherA2.becke_accuracy',beckuracy)
    for basis in basesList:
      OV.SetParam('snum.NoSpherA2.basis_name',basis)
      for epsilon in epsilonList:
        OV.SetParam('lauranonsphertests.vars.epsilon', epsilon)
        olex.m('spy.LauraNonSpherTests.start_data()')      
    
        for j in range(20):
          olex.m('refine')    #refine to aspherical optimum  
          OV.SetParam('snum.NoSpherA2.Calculate',True)       
    
        #Do initial RUN and TEST
        olex.m('spy.LauraNonSpherTests.make_gradient_file()')
        olex.m('spy.LauraNonSpherTests.deposit_phrase("Aspherical Model")')
        olex.m('spy.LauraNonSpherTests.do_test()') 
    
        #Run MULTIPLE
        olex.m('spy.LauraNonSpherTests.deposit_phrase("True Model")')
        olex.m('spy.LauraNonSpherTests.refine_multiple()')


  toc = time.time()
  elapsed_time = toc-tic
  olex.m('spy.LauraNonSpherTests.rename_comparison_file("comparisonData_parameter_choices")')  
  shutil.rmtree(r"olex2\Wfn_job")
  shutil.rmtree(r"olex2\NoSpherA2_history")
  print("Calculation done in "+str(elapsed_time)+" seconds")
  print("PROFIT?!")  

def run_full_test_series(number):
  #does tests for design matrices 1e-1 to 1e-10 for a set number of times
  #not really useful honestly
  import time
  tic = time.time()
  #first: disable NoSpherA2 and refine a few times to make sure things are settled
  OV.SetParam('snum.NoSpherA2.use_aspherical',False)
  olex.m('anis -h')
  olex.m('refine')
  olex.m('refine')
  olex.m('refine')
  olex.m('refine')
  olex.m('refine')

  #put NoSpherA2 settings
  OV.SetParam('snum.NoSpherA2.use_aspherical',True)
  OV.SetParam('snum.NoSpherA2.use_aspherical',True)
  OV.SetParam('snum.NoSpherA2.full_HAR',False)  
  OV.SetParam('snum.NoSpherA2.Calculate',True)  


  #Do initial RUN and TEST
  olex.m('spy.LauraNonSpherTests.make_gradient_file()')
  olex.m('spy.LauraNonSpherTests.start_data()')
  olex.m('spy.LauraNonSpherTests.deposit_phrase("Spherical Model")')
  olex.m('spy.LauraNonSpherTests.do_test()') 


  for i in range(int(number)):

    OV.SetParam('snum.NoSpherA2.use_aspherical',False)
    for i in range(20):
      olex.m('refine')
      olex.m('refine')
      olex.m('refine')
      olex.m('refine')
      olex.m('refine') 

    #Turn on NoSpherA2
    OV.SetParam('snum.NoSpherA2.use_aspherical',True)

    #Run MULTIPLE
    olex.m('spy.LauraNonSpherTests.deposit_phrase("True Model")')
    olex.m('spy.LauraNonSpherTests.refine_multiple()')

    #Go back to spherical model
    OV.SetParam('snum.NoSpherA2.use_aspherical',False)
    olex.m('refine')
    olex.m('refine')
    olex.m('refine')
    olex.m('refine')
    olex.m('refine')

    #Hit RUN and TEST again
    olex.m('spy.LauraNonSpherTests.deposit_phrase("Spherical Model")')
    olex.m('spy.LauraNonSpherTests.make_gradient_file()')
    olex.m('spy.LauraNonSpherTests.do_test()')

  toc = time.time()
  elapsed_time = toc-tic
  print("Calculation done in "+str(elapsed_time)+" seconds")
  print("PROFIT?!")

def run_test_epsilons(number):
  #runs full numerical refinements for 20 steps for each epsilon in the list below (1e-1 to 1e-6 but in order of 'likely usefulness')
  import time
  tic = time.time()
  #first: disable NoSpherA2 and refine a few times to make sure things are settled
  OV.SetParam('snum.NoSpherA2.use_aspherical',False)
  olex.m('anis')
  olex.m('anis -h')
  for i in range(100):
    olex.m('refine')

  #put NoSpherA2 settings
  OV.SetParam('snum.NoSpherA2.use_aspherical',True)
  OV.SetParam('snum.NoSpherA2.full_HAR',False)  
  OV.SetParam('snum.NoSpherA2.source',"DISCAMB")
  OV.SetParam('snum.NoSpherA2.Calculate',True)  
  #OV.SetParam('snum.NoSpherA2.basis_name',"def2-TZVPP")
  OV.SetParam('snum.NoSpherA2.becke_accuracy',"High")
  OV.SetParam('snum.NoSpherA2.method',"PBE")

  olex.m('spy.LauraNonSpherTests.make_gradient_file()')
  olex.m('spy.LauraNonSpherTests.deposit_phrase("Spherical Model")')
  olex.m('spy.LauraNonSpherTests.do_test()') 

  for i in range(20):
    olex.m('refine')    #refine to aspherical optimum 
    OV.SetParam('snum.NoSpherA2.Calculate',True)     

  for i in range(int(number)):
    for epsilon in [1e-4,1e-3,1e-5,1e-2,1e-1,1e-6]:
      OV.SetParam('lauranonsphertests.vars.epsilon', epsilon)
      olex.m('spy.LauraNonSpherTests.start_data()')      

      for j in range(20):
        olex.m('refine')    #refine to aspherical optimum  
        OV.SetParam('snum.NoSpherA2.Calculate',True)       

      #Do initial RUN and TEST
      olex.m('spy.LauraNonSpherTests.make_gradient_file()')
      olex.m('spy.LauraNonSpherTests.deposit_phrase("Aspherical Model")')
      olex.m('spy.LauraNonSpherTests.do_test()') 

      #Run MULTIPLE
      olex.m('spy.LauraNonSpherTests.deposit_phrase("True Model")')
      olex.m('spy.LauraNonSpherTests.refine_multiple()')


  toc = time.time()
  elapsed_time = toc-tic
  olex.m('spy.LauraNonSpherTests.rename_comparison_file("")')  
  shutil.rmtree(r"olex2\Wfn_job")
  shutil.rmtree(r"olex2\NoSpherA2_history")
  print("Calculation done in "+str(elapsed_time)+" seconds")
  print("PROFIT?!")



def run_test_5e4(epsilon=5e-4):
  #runs a test for a particular value of epsilon to record the refienemnt process from spherical, from approximate, and to hybrid.
  import time
  tic = time.time()
  OV.SetParam('lauranonsphertests.vars.epsilon', float(epsilon))
  #first: disable NoSpherA2 and refine a few times to make sure things are settled
  OV.SetParam('snum.NoSpherA2.use_aspherical',False) #this might be toggle-gui now?
  olex.m('anis -h')
  for i in range(20):
    olex.m('refine')

  #put NoSpherA2 settings
  OV.SetParam('snum.NoSpherA2.use_aspherical',True)
  OV.SetParam('snum.NoSpherA2.full_HAR',False)  
  OV.SetParam('snum.NoSpherA2.source',"DISCAMB")
  OV.SetParam('snum.NoSpherA2.Calculate',True)  
  #OV.SetParam('snum.NoSpherA2.basis_name',"def2-TZVPP")
  OV.SetParam('snum.NoSpherA2.becke_accuracy',"High")
  OV.SetParam('snum.NoSpherA2.method',"PBE")

  olex.m('spy.LauraNonSpherTests.start_data()')
  olex.m('spy.LauraNonSpherTests.rename_comparison_file("single_epsilon_parameters")')
  olex.m('spy.LauraNonSpherTests.make_gradient_file()')
  olex.m('spy.LauraNonSpherTests.do_test()') 
  olex.m('spy.LauraNonSpherTests.rename_comparison_file("record_sph")')
  olex.m('spy.LauraNonSpherTests.rename_gradient_file("sph")') 

  olex.m('spy.LauraNonSpherTests.refine_multiple()')  
  olex.m('spy.LauraNonSpherTests.rename_comparison_file("record_tru_from_sph")')

  for i in range(20):
    olex.m('refine')    #refine to aspherical optimum 
    OV.SetParam('snum.NoSpherA2.Calculate',True)     

  olex.m('spy.LauraNonSpherTests.make_gradient_file()')
  olex.m('spy.LauraNonSpherTests.do_test()')   
  olex.m('spy.LauraNonSpherTests.rename_gradient_file("app")')  

  olex.m('spy.LauraNonSpherTests.rename_comparison_file("record_app")')

  #Run MULTIPLE
  olex.m('spy.LauraNonSpherTests.refine_multiple()')
  olex.m('spy.LauraNonSpherTests.rename_gradient_file("tru")')  
  olex.m('spy.LauraNonSpherTests.rename_comparison_file("record_tru_from_app")')


  for i in range(20):
    olex.m('refine')    #refine to aspherical optimum 
    OV.SetParam('snum.NoSpherA2.Calculate',True)      

  olex.m('spy.LauraNonSpherTests.refine_multiple(hybrid=True)')
  olex.m('spy.LauraNonSpherTests.rename_gradient_file("hyb")') 
  olex.m('spy.LauraNonSpherTests.rename_comparison_file("record_hyb_from_app")')

  toc = time.time()
  elapsed_time = toc-tic
  shutil.rmtree(r"olex2\Wfn_job")
  shutil.rmtree(r"olex2\NoSpherA2_history")
  print("Calculation done in "+str(elapsed_time)+" seconds")
  print("PROFIT?!")

def run_all_tests():
  olex.m('spy.LauraNonSpherTests.run_test_epsilons(1)')
  olex.m('spy.LauraNonSpherTests.run_test_5e4(1e-3)')  
  olex.m('spy.LauraNonSpherTests.spherical_refinement_with_ns_fcalc()')

def run_tests_time():
  OV.SetParam('lauranonsphertests.vars.epsilon', 1e-3)
  
  
  nospte=NonSphereTests()
  OV.SetParam('snum.NoSpherA2.Calculate',True)     
  nospte.deposit_phrase("Approx")
  tic=time.time()
  olex.m('refine')  
  toc=time.time()
  nospte = None
  OV.SetParam('snum.NoSpherA2.Calculate',True)     
  elapsed_time=toc-tic
  nospte=NonSphereTests()
  nospte.deposit_phrase(str(elapsed_time))
  nospte.deposit_phrase("Numerical")
  tic=time.time()
  OV.SetParam('snum.NoSpherA2.no_backup',True)
  nospte.run()
  OV.SetParam('snum.NoSpherA2.no_backup',False)
  nospte.refinement_step()
  toc=time.time()
  nospte = None
  elapsed_time=toc-tic
  nospte=NonSphereTests()
  nospte.deposit_phrase(str(elapsed_time))
  nospte.deposit_phrase("Hybrid")
  tic=time.time()
  OV.SetParam('snum.NoSpherA2.no_backup',True)
  nospte.run(hybrid=True)
  OV.SetParam('snum.NoSpherA2.no_backup',False)
  nospte.refinement_step(hybrid=True)
  toc=time.time()
  nospte = None
  elapsed_time=toc-tic
  nospte=NonSphereTests()
  nospte.deposit_phrase(str(elapsed_time))

def many_time():
  for i in range(5):
    olex.m('spy.LauraNonSpherTests.run_tests_time')

OV.registerFunction(run_full_test_series,True,"LauraNonSpherTests")
OV.registerFunction(run_many_tests,True,"LauraNonSpherTests")
OV.registerFunction(run_test_epsilons,True,"LauraNonSpherTests")
OV.registerFunction(run_test_5e4,True,"LauraNonSpherTests")
OV.registerFunction(run_all_tests,True,"LauraNonSpherTests")
OV.registerFunction(run_test_parameter_choices,True,"LauraNonSpherTests")
OV.registerFunction(run_tests_time,True,"LauraNonSpherTests")
OV.registerFunction(many_time,True,"LauraNonSpherTests")