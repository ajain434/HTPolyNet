import logging
import os
# import sys
import argparse as ap
import textwrap

from HTPolyNet.banner import banner, banner_message
from HTPolyNet.runtime import Runtime,logrotate
import HTPolyNet.projectfilesystem as pfs
import HTPolyNet.software as software
from HTPolyNet.plot import cure_graph,density_evolution
from HTPolyNet.stringthings import my_logger

logger=logging.getLogger(__name__)
parser=ap.ArgumentParser()

def info(args):
    print('This is some information on your installed version of HTPolyNet')
    l=pfs.lib_setup()
    software.sw_setup()
    print(l.info())
    print(software.to_string())

def run(args):

    logrotate(args.log)

    ''' set up logger with all debug+-level messages going to the diagnostic log file and info to console '''
    loglevel_numeric=getattr(logging, args.loglevel.upper())
    logging.basicConfig(filename=args.log,filemode='w',format='%(asctime)s %(name)s.%(funcName)s %(levelname)s> %(message)s',level=loglevel_numeric)
    console=logging.StreamHandler()
    console.setLevel(logging.INFO)
    formatter=logging.Formatter('%(levelname)s> %(message)s')
    console.setFormatter(formatter)
    logging.getLogger('').addHandler(console)

    banner(logger.info)
    my_logger('HTPolyNet runtime begins',logger.info)
    userlib=args.lib
    if not os.path.exists(args.lib):
        userlib=None
    software.sw_setup()
    pfs.pfs_setup(root=os.getcwd(),topdirs=['molecules','systems','plots'],verbose=True,projdir=args.proj,reProject=args.restart,userlibrary=userlib)
    a=Runtime(cfgfile=args.config,restart=args.restart)
    a.build(force_checkin=args.force_checkin,force_parameterization=args.force_parameterization)
    my_logger('HTPolyNet runtime ends',logger.info)

def parameterize(args):

    logrotate(args.log)

    ''' set up logger with all debug+-level messages going to the diagnostic log file and info to console '''
    loglevel_numeric=getattr(logging, args.loglevel.upper())
    logging.basicConfig(filename=args.log,filemode='w',format='%(asctime)s %(name)s.%(funcName)s %(levelname)s> %(message)s',level=loglevel_numeric)
    console=logging.StreamHandler()
    console.setLevel(logging.INFO)
    formatter=logging.Formatter('%(levelname)s> %(message)s')
    console.setFormatter(formatter)
    logging.getLogger('').addHandler(console)

    banner(logger.info)
    my_logger('HTPolyNet parameterization begins',logger.info)
    userlib=args.lib
    if not os.path.exists(args.lib):
        userlib=None
    software.sw_setup()
    pfs.pfs_setup(root=os.getcwd(),topdirs=['molecules','systems','plots'],verbose=True,reProject=args.restart,userlibrary=userlib)
    a=Runtime(cfgfile=args.config,restart=args.restart)
    a.generate_molecules(force_checkin=args.force_checkin,force_parameterization=args.force_parameterization)
    my_logger('HTPolynet parameterization ends',logger.info)

def htpolynet_cure_plots(args):
    logs=args.logs
    banner(print)
    cure_graph(logs,args.plotfile)
    density_evolution()

def cli():
    """cli Command-line interface
    """
    commands={}
    commands['run']=run
    commands['parameterize']=parameterize
    commands['info']=info
    commands['plots']=htpolynet_cure_plots

    helps={}
    helps['run']='build a system using instructions in the config file and any required molecular structure inputs'
    helps['parameterize']='parameterize monomers and oligomer templates using instructinos in the config file'
    helps['info']='print some information to the console'
    helps['plots']='generate some plots that summarize aspects of the current completed build'

    parser=ap.ArgumentParser(description=textwrap.dedent(banner_message),formatter_class=ap.RawDescriptionHelpFormatter)
    subparsers=parser.add_subparsers()
    command_parsers={}
    for k in commands:
        command_parsers[k]=subparsers.add_parser(k,help=helps[k])
        command_parsers[k].set_defaults(func=commands[k])

    command_parsers['run'].add_argument('config',type=str,default=None,help='input configuration file in YAML format')
    command_parsers['run'].add_argument('-lib',type=str,default='lib',help='local user library of molecular structures and parameterizations')
    command_parsers['run'].add_argument('-proj',type=str,default='next',help='project directory; "next" (default) generates next directory\nAnything other than "next": if it exists, "-restart" must be included as a parameter; if not, it is created as a new project')
    command_parsers['run'].add_argument('-log',type=str,default='htpolynet_runtime_diagnostics.log',help='diagnostic log file')
    command_parsers['run'].add_argument('-restart',default=False,action='store_true',help='restart in latest proj dir')
    command_parsers['run'].add_argument('--force-parameterization',default=False,action='store_true',help='force GAFF parameterization of any input mol2 structures')
    command_parsers['run'].add_argument('--force-checkin',default=False,action='store_true',help='force check-in of any generated parameter files to the system library')
    command_parsers['run'].add_argument('--loglevel',type=str,default='debug',help='Log level for messages written to diagnostic log (debug|info)')

    command_parsers['parameterize'].add_argument('config',type=str,default=None,help='input configuration file in YAML format')
    command_parsers['parameterize'].add_argument('-lib',type=str,default='lib',help='local user library of molecular structures and parameterizations')
    command_parsers['parameterize'].add_argument('-log',type=str,default='htpolynet_runtime_diagnostics.log',help='diagnostic log file')
    command_parsers['parameterize'].add_argument('-restart',default=False,action='store_true',help='restart in latest proj dir')
    command_parsers['parameterize'].add_argument('--force-parameterization',default=False,action='store_true',help='force GAFF parameterization of any input mol2 structures')
    command_parsers['parameterize'].add_argument('--force-checkin',default=False,action='store_true',help='force check-in of any generated parameter files to the system library')
    command_parsers['parameterize'].add_argument('--loglevel',type=str,default='debug',help='Log level for messages written to diagnostic log (debug|info)')

    command_parsers['plots'].add_argument('logs',type=str,default='',nargs='+',help='names of diagnostic log files')
    command_parsers['plots'].add_argument('--plotfile',type=str,default='cure-info.png',help='name of plot file to generate')

    args=parser.parse_args()
    args.func(args)
