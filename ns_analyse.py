import numpy as np
import copy

#startup details
cellParamEpox=np.array([4.633, 8.4, 6.577, 90, 100.37, 90])
cellSigmaEpox=np.array([0.005,0.001, 0.003, 0, 0.06, 0])
cellParamAlan=np.array([5.791,5.944,12.269,90,90,90])
cellSigmaAlan=np.array([0.002,0.002, 0.002, 0, 0, 0])
cellParamAmm=np.array([5.1305,5.1305,5.1305,90,90,90]) 
cellSigmaAmm=np.array([0.0008,0.0008, 0.0008, 0, 0, 0])
cellParamAlan[3:]*=np.pi/180
cellSigmaAlan[3:]*=np.pi/180
cellParamAmm[3:] *=np.pi/180
cellSigmaAmm[3:] *=np.pi/180
cellParamEpox[3:]*=np.pi/180
cellSigmaEpox[3:]*=np.pi/180
orthogEpox=np.array([[ 4.633,  5.14351656e-16, -1.18388712], [ 0.,  8.4,  1.04038626e-15],  [ 0.,  0.,  6.46957033]])
orthogAla=np.array([[5.791, 3.63965029e-16, 7.51259579e-16],  [0., 5.944, 1.97300568e-15],  [0., 0., 1.2269e+01]])
orthogAmm=np.array([[5.1305,3.14152520e-16,3.14152520e-16],  [0, 5.1305, 8.25047325e-16],  [0, 0, 5.1305]])
hklEpox="epoxide2.hkl"
hklAmm="nh30r2.hkl"
hklAla="alanine2.hkl"


class record_set:
    def __init__(self, appModel=None, truModel=None, inputDetails=None,steps=None):
        #input details, approx model, true model
        if not steps:
            steps=[]
        self.appModel=appModel
        self.truModel=truModel
        if inputDetails:
            self.source=inputDetails[0]
            self.calculate=inputDetails[1]
            self.basis=inputDetails[2]
            self.becke=inputDetails[3]
            self.method=inputDetails[4]
            self.har=inputDetails[5]
            self.eps=inputDetails[6]
        self.steps=steps
    def addStep(self, step):
        self.steps+=[step]
        
class structure:
    def __init__(self, atoms=None, varCovar=None):
        #list of atoms and locations, var/covar
        if not atoms:
            atoms=[]
        if not varCovar:
            varCovar=[]
        self.atoms=atoms
        self.varCovar=varCovar
        self.gradient_records=None
    def addAtom(self, atom):
        self.atoms+=[atom]
    def addCovarLine(self,line):
        self.varCovar+=[line]
    def setRfac(self,rfacs):
        self.r1cut   , self.wr2cut  = np.fromstring((rfacs.rsplit("(")[1]).rsplit(")")[0],sep=',')
        self.r1nocut , self.wr2nocut= np.fromstring((rfacs.rsplit("(")[2]).rsplit(")")[0],sep=',')
    def getParameterList(self):
        paramList=[]
        for atom in self.atoms:
            paramList+=[atom.params]
        return paramList
    def getPosList(self):
        paramList=[]
        for atom in self.atoms:
            paramList+=[atom.pos]
        return paramList
    def getAdpList(self):
        paramList=[]
        for atom in self.atoms:
            paramList+=[atom.adp]
        return paramList
    def structure_difference(self,struct2):
        posList=[]
        adpList=[]
        for i, atom in enumerate(self.atoms):
            atom2=struct2.atoms[i]
            euclidLine=atom.params_as_euclid()-atom2.params_as_euclid()
            posList+=list(euclidLine[:3])
            adpList+=list(euclidLine[3:])
        return(posList,adpList)
    def add_gradient_records(self,filename):
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
        orthog=[]
        params=[]
    
        lines = test.readlines()
    
        i=0
        for line in lines:
            if line[0]!='[' and line[0]!=' ':
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
            elif i==21:
                arrayline=np.fromstring((line.rsplit("]")[0]).rsplit("[")[-1],sep=' ')
                orthog+=[list(arrayline)]
            elif i==22:
                params+=[line]
    
        test.close()        
    
        grads_aspher_tru=np.array(grads_aspher_tru)
        grads_aspher_app=np.array(grads_aspher_app)
        grads_aspher_hybrid=np.array(grads_aspher_hybrid)
        grads_spher=np.array(grads_spher)
        tsc_design_ma=np.array(tsc_design_ma)
        design_ma=np.array(design_ma) 
        orthog=np.array(orthog)
        fcalcs_spher=fcalcs_spher[0]
        fcalcs_tru=fcalcs_tru[0]
        weights_spher=weights_spher[0]
        weights_tru=weights_tru[0]
    
        self.gradient_records=(grads_aspher_tru,grads_aspher_app,grads_aspher_hybrid,grads_spher,tsc_design_ma,design_ma,fcalcs_tru,weights_tru,fcalcs_spher,weights_spher,orthog,params)    
        
    
class step:
    def __init__(self, line1, line2):
        #max pos shift, max adp shift, r factors
        self.pos     , self.adp     = np.fromstring(line1.rsplit("=")[1],sep=',')
        self.r1cut   , self.wr2cut  = np.fromstring((line2.rsplit("(")[1]).rsplit(")")[0],sep=',')
        self.r1nocut , self.wr2nocut= np.fromstring((line2.rsplit("(")[2]).rsplit(")")[0],sep=',')
    def getRow(self):
        row=str(self.pos)+" "+str(self.adp)+" "+str(self.r1cut*100)+" "+str(self.wr2nocut*100-1.2)#-2.68)#-5.26452)# #BOTCH ZZZZZ )#)
        return row

class atom:
    def __init__(self, atomDataLine):
        parts=atomDataLine.rsplit("(")
        self.atomName=parts[0][2:-3]
        self.pos=np.fromstring((parts[1].rsplit(")"))[0],sep=',')
        self.adp=np.fromstring((parts[2].rsplit(")"))[0],sep=',')
        self.params=np.concatenate((self.pos,self.adp))
    def params_as_euclid(self):
        return euclid_line(self.params)
    def params_as_olex_display(self):
        return olex_display(self.params)
    def distance_from(self, atom2):
        difference=self.params[0:3]-atom2.params[0:3]
        euclid_difference=np.dot(orthog,difference)
        size=np.linalg.norm(euclid_difference)
        return size
    def distance_from_adp(self, atom2):
        difference=self.params[3:]-atom2.params[3:]
        euclid_difference=u_rel_to_ucart(difference)
        size=np.linalg.norm(euclid_difference)
        return size

class Molecule_set:
    def __init__(self,name,orthog,cellParam,cellSigma,n_atoms,epsChoice=None):
        self.molName=name
        self.orthog=orthog
        self.orthogInv=np.linalg.inv(orthog)
        self.metric=np.dot(orthog.T,orthog)   
        self.cellParam=cellParam
        self.cellSigma=cellSigma
        self.n_atoms=n_atoms
        self.manyEpsilonRecordsList=None
        self.spherBase=None
        self.truFromSpher=None
        self.truFromApp=None
        self.mixFromSpher=None
        self.hybFromApp=None
        self.sphRecord=None
        self.appRecord=None
        self.epsChoice=None
    
    def retrieve_single_record(self,filename):
        record=record_set()
        file=open(filename)
        lines=file.readlines()
        n=self.n_atoms
        k=self.n_atoms*10
        for i in range(21):
            record.addStep(step(lines[2*i], lines[2*i+1]))
        theStructure=structure()
        for i in range(43,43+n):
            theStructure.addAtom(atom(lines[i]))
        for i in range(43+n,43+k):
            theStructure.addCovarLine(np.fromstring(lines[i][1:-1],sep=","))
        theStructure.setRfac(lines[43+k])
        record.truModel=theStructure
        return record    

    def retrieve_single_limited_record(self,filename):
        record=record_set()
        file=open(filename)
        lines=file.readlines()
        theStructure=structure()
        n=self.n_atoms
        k=self.n_atoms*10
        for i in range(0,n):
            theStructure.addAtom(atom(lines[i]))
        for i in range(n,k):
            theStructure.addCovarLine(np.fromstring(lines[i][1:-1],sep=","))
        theStructure.setRfac(lines[k])
        record.truModel=theStructure
        return record    
    
    
    def retrieve_records_from_full_records(self,filename):
        recordsList=[]
        spherBase=structure()
        
        currentStructure=None
        currentRecord=None
        
        file=open(filename)
        lines=file.readlines()
        total=len(lines)
        i=0
        done=False
        for j in range(len(lines)): #using a for so that we can be sure it will definitely exit if I mess up.
            if done==True:
                break
            iIncreased=False
            line=lines[i]
            if line[0]=="[":
                #this is the case when we are on a structure line X, or the uncertainty line X, or the r-factors X. It can also be the start of the orthogonalisation matrix.
                if line[1]=="\'":
                    #we are in the model case or the R factor case
                    if line[2]=="R":
                        currentStructure.setRfac(line)
                    else:
                        while line[1]=="\'":
                            newAtom=atom(line)
                            currentStructure.addAtom(newAtom)
                            i+=1
                            if i<total:
                                line=lines[i]
                                iIncreased=True
                            else:
                                done=True
                                break
                elif line[1]=="[":
                    orthog=[]
                    orthog+=[np.fromstring(line.rsplit("[")[-1][:-1],sep=" ")]
                    i=i+1
                    line=lines[i]
                    orthog+=[np.fromstring(line.rsplit("[")[-1][:-1],sep=" ")]
                    i=i+1
                    line=lines[i]                
                    orthog+=[np.fromstring(line.rsplit("]")[0][2:],sep=" ")]
                else:
                    #we are in the uncertainty line
                    while line[1]!="\'":
                        currentStructure.addCovarLine(line)
                        i=i+1
                        if i<total:
                            line=lines[i]
                            iIncreased=True
                        else:
                            done=True
                            break
            elif line=="New Record\n":
                if currentRecord:
                    currentRecord.truModel=copy.deepcopy(currentStructure)
                    recordsList+=[copy.deepcopy(currentRecord)]
                    currentRecord=None
                    currentStructure=None
                else:
                    spherBase=copy.deepcopy(currentStructure)
                    currentStructure=None
                currentRecord=record_set()
                i+=2
                if i<total:
                    line=lines[i]
                    iIncreased=True
                else:
                    done=True
                    break                
            elif line[:3]=="Sou":
                currentRecord.source=line[8:]
            elif line[:3]=="Cal":
                currentRecord.calculate=line[11:]   
            elif line[:3]=="Bas":
                currentRecord.basis=line[12:]
            elif line[:3]=="Bec":
                currentRecord.becke=line[16:]
            elif line[:3]=="Met":
                currentRecord.method=line[8:]
            elif line[:3]=="Ful":
                currentRecord.har=line[10:]
            elif line[:3]=="Eps":
                currentRecord.eps=float(line[9:-1])
            elif line[-6:-1]=="Model":
                if line[:3]=="Sph":
                    currentStructure=structure()
                elif line[:3]=="Asp":
                    currentStructure=structure()
                elif line[:3]=="Tru":
                    currentRecord.appModel=copy.deepcopy(currentStructure)
                    currentStructure=None
                    currentStructure=structure()
                    i=i+1
                    if i<total:
                        line=lines[i] 
                    else:
                        done=True
                        break
                    while line[:5]!="Steps":
                        currentRecord.addStep(step(line,lines[i+1]))
                        i=i+2
                        if i<total:
                            line=lines[i]
                            iIncreased=True
                        else:
                            done=True
                            break
                    i=i+1
                    if i<total:
                        line=lines[i] 
                        iIncreased=True
                    else:
                        done=True
                        break
            else:
                print("Line Not Found: "+line)
            if not iIncreased:
                i=i+1
                if i<total:
                    line=lines[i] 
                else:
                    done=True
                    break
            #finally
            
        self.manyEpsilonRecordsList=recordsList
        self.spherBase=spherBase
        #return (recordsList,spherBase)              
    
    def populate_consistent_eps_test(self):
        name=self.molName
        if self.epsChoice:
            name+=self.epsChoice
        self.truFromSpher=self.retrieve_single_record("record_tru_from_spher_"+name+".txt")
        self.truFromApp=self.retrieve_single_record("record_tru_from_app_"+name+".txt")
        self.hybFromApp=self.retrieve_single_record("record_hyb_from_app_"+name+".txt") 
        self.sphRecord=self.retrieve_single_limited_record("record_sph_"+name+".txt")
        self.appRecord=self.retrieve_single_limited_record("record_app_"+name+".txt")
        self.mixFromSpher=self.retrieve_single_record("record_mix_from_sph_"+name+".txt")
        
    def populate_epsilon_tests(self):
        name=self.molName
        self.retrieve_records_from_full_records("comparisonData_"+name+".txt")
        

def create_data():
    ammBase=Molecule_set("nh30r", orthogAmm, cellParamAmm, cellSigmaAmm, 2)
    alaBase=Molecule_set("alanine", orthogAla, cellParamAlan, cellSigmaAlan, 13)
    epoBase=Molecule_set("epoxide", orthogEpox, cellParamEpox, cellSigmaEpox, 7)