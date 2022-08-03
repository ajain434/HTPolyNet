import os
import yaml
import pandas as pd
from enum import Enum
import logging
from HTPolyNet.topocoord import TopoCoord

logger=logging.getLogger(__name__)

class state(Enum):
    """Enumerated CURE state
    """
    bondsearch=0
    drag=1
    update=2
    relax=3
    equilibrate=4
    postcure=5
    postcure_equilibrate=6
    finished=7
    unknown=99
    def __str__(self):
        return self.name

class Checkpoint:
    def __init__(self):
        self.state=state.unknown
        self.iter=0
        self.curr_conversion=0.0
        self.max_search_radius=0.0
        self.max_radidx=0
        self.curr_nxlinkbonds=0
        self.max_nxlinkbonds=0
        self.current_dragstage=0
        self.current_stage=0
        self.current_radidx=0
        self.radius=0.0
        self.checkpoint_file=''
        self.bonds_file=None
        self.bonds=pd.DataFrame()
        self.bonds_are='nonexistent!'
        self.top=None
        self.gro=None
        self.grx=None
        self.mol2=None
        self.cwd=os.getcwd()
        self.stepschecked=[]

    def pfx(self):
        return f'{self.state.value}-{self.state}'

    def set_state(self,state):
        self.state=state

    def is_cured(self):
        return self.curr_nxlinkbonds==self.max_nxlinkbonds

    def reset_for_next_cure_iter(self):
        self.iter+=1
        self.state=state.bondsearch
        self.curr_conversion=0.0
        self.max_search_radius=0.0
        self.max_radidx=0
        self.curr_nxlinkbonds=0
        self.current_dragstage=0
        self.current_stage=0
        self.current_radidx=0
        self.radius=0.0
        self.bonds_file=None
        self.bonds=pd.DataFrame()
        self.bonds_are='nonexistent!'

    def _write_bondsfile(self):
        self.bonds.to_csv(self.bonds_file,sep=' ',mode='w',index=False,header=True,doublequote=False)

    def _read_bondsfile(self):
        assert os.path.exists(self.bonds_file),f'Error: {self.bonds_file} not found.'
        self.bonds=pd.read_csv(self.bonds_file,sep='\s+',header=0)

    def register_bonds(self,bonds,pairs,bonds_are='unrelaxed'):
        self.bonds=bonds
        self.pairs=pairs
        self.bonds_are=bonds_are
    
    def read_checkpoint(self): #,system):
        if os.path.exists(self.checkpoint_file):
            with open(self.checkpoint_file,'r') as f:
                basedict=yaml.safe_load(f)
            self.top=os.path.basename(basedict['topology'])
            self.gro=os.path.basename(basedict['coordinates'])
            self.grx=os.path.basename(basedict['extra_attributes'])
            self.mol2=basedict.get('mol2_coordinates',None)
            if self.mol2:
                self.mol2=os.path.basename(basedict['mol2_coordinates'])
            self.cwd=os.path.commonpath([basedict['topology'],basedict['coordinates']])
            os.chdir(self.cwd)
            # self.iter=basedict['ITERATION']
            # self.state=CPstate[basedict['STATE']]
            # self.current_dragstage=basedict['CURRENT_DRAGSTAGE']
            # self.current_stage=basedict['CURRENT_STAGE']
            # self.current_radidx=basedict['CURRENT_RADIDX']
            # self.radius=basedict['RADIUS']
            # bf=basedict.get('BONDSFILE',None)
            # # assert bf,f'Error: BONDSFILE not found in {self.checkpoint_file}.'
            # self.bonds_are=basedict.get('BONDS_ARE',None)
            # if bf:
            #     self.bonds_file=os.path.basename(bf)
            #     self._read_bondsfile()
            # else:
                # self.bonds_file=None
            #system.set_system(CP=self)
        # else:
        #     logging.debug(f'read_checkpoint: no file, empty checkpoint')

    def write_checkpoint(self,TC:TopoCoord,stepname): #,state): #system,state,prefix='checkpoint'):
        # self.state=state
        self.stepschecked.append(stepname)
        self.top,self.gro,self.grx,self.mol2=[os.path.basename(TC.files[x]) for x in ['top','gro','grx','mol2']]
        if self.top:
            TC.write_top(self.top)
        if self.gro:
            TC.write_gro(self.gro)
        if self.grx:
            TC.write_grx_attributes(self.grx)
        if self.mol2:
            TC.write_mol2(self.mol2)
        self.cwd=os.getcwd()
        # self.bonds_file=prefix+'-bonds.csv'
        # system.register_system(CP=self)
        with open(self.checkpoint_file,'w') as f:
            # f.write(f'ITERATION: {self.iter}\n')
            # f.write(f'STATE: {str(self.state)}\n')
            # f.write(f'CURRENT_DRAGSTAGE: {self.current_dragstage}\n')
            # f.write(f'CURRENT_STAGE: {self.current_stage}\n')
            # f.write(f'CURRENT_RADIDX: {self.current_radidx}\n')
            # f.write(f'RADIUS: {self.radius}\n')
            f.write(f'stepschecked: {self.stepschecked}\n')
            f.write(f'cwd: {self.cwd}\n')
            if self.top:
                f.write(f'topology: {self.top}\n')
            if self.gro:
                f.write(f'coordinates: {self.gro}\n')
            if self.grx:
                f.write(f'extra_attributes: {self.grx}\n')
            if self.mol2:
                f.write(f'mol2_coordinates: {self.mol2}\n')
            # if self.bonds.shape[0]>0:
            #     f.write(f'BONDS_ARE: {self.bonds_are}\n')
            #     f.write(f'BONDSFILE: {self.bonds_file}\n')
            f.close()
        # self._write_bondsfile()

_CP_=Checkpoint()
def setup(checkpoint_file='checkpoint.yaml'):
    global _CP_
    if os.path.exists(checkpoint_file):
        _CP_.checkpoint_file=os.path.abspath(checkpoint_file)
    else:
        _CP_.checkpoint_file=os.path.join(os.getcwd(),checkpoint_file)
    logger.debug(f'checkpoints to {_CP_.checkpoint_file}')

def write(TC:TopoCoord,stepname):#,state:CPstate):
    global _CP_
    _CP_.write_checkpoint(TC,stepname)

def read():
    global _CP_
    _CP_.read_checkpoint()
    return TopoCoord(topfilename=_CP_.top,grofilename=_CP_.gro,grxfilename=_CP_.grx,mol2filename=_CP_.mol2)

def passed(stepname):
    return stepname in _CP_.stepschecked
