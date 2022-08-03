''' Check for presence of required software '''
import subprocess
import logging
import os
logger=logging.getLogger(__name__)

class Software:
    ambertools=['antechamber','tleap','parmchk2']
    def __init__(self):
        cnf=[]
        passes=True
        for c in Software.ambertools:
            CP=subprocess.run(['which',c],capture_output=True,text=True)
            if CP.returncode!=0:
                passes=False
                cnf.append(c)
        assert passes,f'Could not find {cnf}'

    def set_gmx_preferences(self,parameters):
        gromacs_dict=parameters.get('gromacs',{})
        logger.debug(f'gromacs_dict {gromacs_dict}')
        if gromacs_dict:
            self.gmx=gromacs_dict.get('gmx','gmx')
            self.gmx_options=gromacs_dict.get('gmx_options','-quiet')
            self.mdrun=gromacs_dict.get('mdrun','gmx mdrun')
            logger.debug(f'{self.gmx}, {self.gmx_options}, {self.mdrun}')
        else:
            self.gmx_options=parameters.get('gmx_options','')
            self.gmx=parameters.get('gmx','gmx')
            self.mdrun=parameters.get('gmx_mdrun',f'{self.gmx} {self.gmx_options} mdrun')
        CP=subprocess.run(['which',self.gmx],capture_output=True,text=True)
        assert CP.returncode==0,f'{self.gmx} not found'

    def __str__(self):
        self.getVersions()
        r='Ambertools commands available for HTPolyNet to use:\n'
        for c in self.ambertools:
            r+=f'{os.path.split(c)[1]:>12s} (ver. {self.versions["ambertools"]:>6s}) at {c:<50s}\n'
        return r

    def getVersions(self):
        self.versions={}
        CP=subprocess.run(['antechamber','-h'],capture_output=True,text=True)
        l=CP.stdout.split('\n')[1].split()[3].strip().strip(':')
        self.versions['ambertools']=l
        
    def info(self):
        print(str(self))

_SW_:Software=None
def sw_setup():
    global _SW_
    _SW_=Software()

def info():
    _SW_.info()

def to_string():
    return str(_SW_)

gmx='gmx'
gmx_options='-quiet'
mdrun=f'{gmx} mdrun'
def set_gmx_preferences(parmdict):
    global gmx
    global gmx_options
    global mdrun
    global _SW_
    _SW_.set_gmx_preferences(parmdict)
    gmx=_SW_.gmx
    gmx_options=_SW_.gmx_options
    mdrun=_SW_.mdrun