"""
openMM code to calculate osmotic pressure and osmotic coefficients from harmonic potentials

"""

# Imports
import logging
logging.basicConfig(
    level=logging.INFO
)

import warnings
warnings.filterwarnings('ignore')

from typing import Union, Iterable, Optional

import argparse
from argparse import Namespace

import os
import os.path
import shutil
import pickle
from pathlib import Path
import numpy as np
import pandas as pd
import math
import json
import csv
from dataclasses import dataclass
from tqdm import tqdm
import subprocess
import shutil
from builtins import sum
from tqdm.notebook import tqdm 

import openmm
from openmm import CustomExternalForce, System, CustomCentroidBondForce, CustomNonbondedForce
from openmm.app import Topology as OMMTopology
from openmm import MonteCarloAnisotropicBarostat,Vec3

from openff.toolkit import Molecule, Topology
from openff.toolkit import ForceField as ForceField
from openff.interchange import Interchange

from polymerist.genutils.fileutils.pathutils import assemble_path
from polymerist.genutils.decorators.functional import allow_string_paths

from polymerist.mdtools.openfftools import topology
from polymerist.mdtools.openfftools import boxvectors
from polymerist.mdtools.openfftools import TKWRAPPERS, GTR

from polymerist.mdtools.openfftools.unitsys import openff_to_openmm
from polymerist.mdtools.openfftools.solvation.solvents import water_TIP3P
from polymerist.mdtools.openfftools.partialcharge.molchargers import NAGLCharger
from polymerist.mdtools.openfftools import FF_PATH_REGISTRY, FF_DIR_REGISTRY

from polymerist.mdtools.openmmtools.execution import run_simulation_schedule
from polymerist.mdtools.openmmtools.parameters import SimulationParameters

from openff.units import UnitRegistry
ureg = UnitRegistry()

from polymerist.mdtools.openfftools.partialcharge.molchargers import NAGLCharger, EspalomaCharger, ABE10Charger
from openff.toolkit.utils.nagl_wrapper import NAGLToolkitWrapper
from openff.toolkit.utils.openeye_wrapper import OpenEyeToolkitWrapper
from openff.toolkit.utils.rdkit_wrapper import RDKitToolkitWrapper
from openff.toolkit.utils.ambertools_wrapper import AmberToolsToolkitWrapper
from openff.toolkit.utils.builtin_wrapper import BuiltInToolkitWrapper

import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.ticker import PercentFormatter

from scipy.optimize import least_squares, minimize
from scipy.integrate import simpson
from scipy.integrate import quad, simpson, trapezoid
from timeit import default_timer as timer

# from openmm.unit import *
from openmm.unit import bar, mole, litre, kelvin, kilojoule_per_mole, nanometer, angstrom, kilocalorie_per_mole, kilogram, molar, atmosphere, nanosecond, picosecond, femtoseconds
from openmm.unit import Quantity, Unit
from openmm.unit import AVOGADRO_CONSTANT_NA, BOLTZMANN_CONSTANT_kB

import MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array
from MDAnalysis import units
from MDAnalysis.analysis import rdf
from MDAnalysis.analysis import density
from pymbar import timeseries

from polymerist.genutils.fileutils.pathutils import assemble_path
from polymerist.mdtools.openmmtools.parameters import SimulationParameters, ThermoParameters, IntegratorParameters, ReporterParameters
from polymerist.mdtools.openmmtools.thermo import EnsembleFactory

from salt_data import SaltData, salt_infos


## Loading information for specified salt
def load_salt_info(ion1, ion2):
    """
    Load all data entries for a specific salt from salt_info and return them as a dictionary
    where the salt appears only once with a subdictionary containing all its data.
    """
    salt = ion1 + ion2
    print(f"Salt to be analyzed: {salt}")

    # Convert dictionaries into SaltData instances
    salt_info_cleaned = [SaltData(**entry) if isinstance(entry, dict) else entry for entry in salt_infos]

    # Filter all entries that match the requested salt and store in a subdictionary
    filtered_entries = {
        f"Molality {entry.molality} mol/kg": {
            "Molarity": entry.molarity,
            "Number of Particles": math.ceil(entry.num_particles),
            "Osmotic Coefficient": entry.osmotic_coefficient,
        }
        for entry in salt_info_cleaned if entry.salt == salt
    }

    # Return results with the salt name as the top-level key
    if filtered_entries:
        return {salt: filtered_entries}
    else:
        return {"Error": f"No data found for {salt}"}
    
def modify_omm_flat_bottom(concentration,repnum,input_picklefile,ion1,ion2,wdir,center_atom1,center_atom2,k,delta_z,z_center):
    with open(input_picklefile, "rb") as f:
        omm_build = pickle.load(f)

    omm_objects_mod={}
    
    for r in tqdm(range(repnum), desc=f"Modifying {concentration}m systems", leave=False):
        rep_key = f'r{r}'
        print(rep_key)

        omm_top = omm_build[rep_key]["topology"]
        omm_pos = omm_build[rep_key]["positions"]
        omm_sys = omm_build[rep_key]["system"]

        # Step 1: Identify polyatomic center atoms for both ions
        poly_atoms = []

        for residue in omm_top.residues():
            if residue.name == ion1.upper():  
                for atom in residue.atoms():
                    if atom.element.symbol == f"{center_atom1}":  # Find the phosphorus atom
                        poly_atoms.append(atom.index)
            if residue.name == ion2.upper():  
                for atom in residue.atoms():
                    if atom.element.symbol == f"{center_atom2}":  # Find the phosphorus atom
                        poly_atoms.append(atom.index)

        print("CENTER ATOMS",len(poly_atoms),poly_atoms)         

        # Step 2: Apply CustomExternalForce to dummy particles
        fb_force = CustomExternalForce('0.5*k*(max(0, abs(z-z0)-rbf)^2)')
        fb_force.addGlobalParameter('k', k)
        fb_force.addGlobalParameter('rbf', delta_z)
        fb_force.addGlobalParameter('z0', z_center)

        for atom_index in poly_atoms:
            fb_force.addParticle(atom_index, [])

        omm_sys.addForce(fb_force)

        # Step 3: Add anisotropic barostat to see if higher concentrations are possible for polyatomics
        pressure_vector = Vec3(1.01325, 1.01325, 0.0) * atmosphere
        anisotropic_force=MonteCarloAnisotropicBarostat(pressure_vector,300,True,True,False)

        omm_sys.addForce(anisotropic_force)

        # Save each replicate's results in the dictionary
        omm_objects_mod[f"r{r}"] = {
            'topology': omm_top,
            'positions': omm_pos,
            'system': omm_sys
        }
        
    # Save the full dictionary
    output_picklefile = f"{wdir}/omm_modified_{concentration}.pkl"
    with open(output_picklefile, "wb") as f:
        pickle.dump(omm_objects_mod, f)
    print(f"✅ System dictionary for {concentration}m saved as pickle file")

    return omm_objects_mod
def build_system_FBP(concentration, repnum, wdir, ion1, ion2, ff, water, salt_dict,center_atom1,center_atom2,k,delta_z,z_center):
    
    sdf_path1=f'structures/{ion1.lower()}.sdf'
    sdf_path2=f'structures/{ion2.lower()}.sdf'

    if concentration % 1 == 0:
        mi1 = f"{concentration:.1f}"
        mi = str(int(concentration))  # Whole number: strip decimal
    else:
        mi = f"{concentration:.1f}".replace('.', '')  # Float: remove the dot
        mi1 = f"{concentration:.1f}"  # used in salt_dict key lookup

    print(f'Molality of maximum concentration = {mi} mol/kg',
            salt_dict[f'{ion1}{ion2}'][f'Molality {mi1} mol/kg'])
    
    input_picklefile=f"{wdir}/omm_build_{mi}.pkl"

    if os.path.exists(input_picklefile):
        print("File already exists. Skipping code.")
        with open(input_picklefile, "rb") as f:
            omm_build = pickle.load(f)
        return modify_omm_flat_bottom(mi,repnum,input_picklefile,ion1,ion2,wdir,center_atom1,center_atom2,k,delta_z,z_center)
   
    else:
        print("File not found. Running code...")

        omm_builds = {}

        for r in tqdm(range(repnum), desc=f"Building {mi}m systems", leave=False):
            pdb_path = f'structures/{ion1.lower()}{ion2.lower()}_{mi}m_r{r}.pdb'
            POL1 = Molecule.from_file(sdf_path1)
            POL2 = Molecule.from_file(sdf_path2)
            
            off_top = Topology.from_pdb(pdb_path, unique_molecules=[POL1,POL2])

            if water == 'TIP3P':
                inc = ff.create_interchange(
                    topology=off_top,
                    toolkit_registry=GTR,
                    charge_from_molecules=[water_TIP3P]
                )
            else:
                inc = ff.create_interchange(
                    topology=off_top,
                    toolkit_registry=GTR
                )

            inc.box = boxvectors.get_topology_bbox(off_top)

            omm_top = inc.to_openmm_topology(collate=True)
            omm_pos = openff_to_openmm(inc.get_positions(include_virtual_sites=True))
            omm_sys = inc.to_openmm_system(combine_nonbonded_forces=False, add_constrained_forces=True)

            # Save each replicate's results in the dictionary
            omm_builds[f"r{r}"] = {
                'topology': omm_top,
                'positions': omm_pos,
                'system': omm_sys
            }

        # Save the full dictionary
        
        with open(input_picklefile, "wb") as f:
            pickle.dump(omm_builds, f)
        print(f"✅ System dictionary for {mi1}m saved as pickle file")
        
        return modify_omm_flat_bottom(mi,repnum,input_picklefile,ion1,ion2,wdir,center_atom1,center_atom2,k,delta_z,z_center)
    

def modify_omm_harmonic(concentration,repnum,input_picklefile,ion1,ion2,wdir,center_atom1,center_atom2,k,delta_z,z_center):
    with open(input_picklefile, "rb") as f:
        omm_build = pickle.load(f)

    omm_objects_mod={}
    
    for r in tqdm(range(repnum), desc=f"Modifying {concentration}m systems", leave=False):
        rep_key = f'r{r}'
        print(rep_key)

        omm_top = omm_build[rep_key]["topology"]
        omm_pos = omm_build[rep_key]["positions"]
        omm_sys = omm_build[rep_key]["system"]

        # Step 1: Identify polyatomic center atoms for both ions
        poly_atoms = []

        for residue in omm_top.residues():
            if residue.name == ion1.upper():  
                for atom in residue.atoms():
                    if atom.element.symbol == f"{center_atom1}":  # Find the phosphorus atom
                        poly_atoms.append(atom.index)
            if residue.name == ion2.upper():  
                for atom in residue.atoms():
                    if atom.element.symbol == f"{center_atom2}":  # Find the phosphorus atom
                        poly_atoms.append(atom.index)

        print("CENTER ATOMS",len(poly_atoms),poly_atoms)    

        # Step 2: Apply CustomExternalForce to dummy particles
        fb_force = CustomExternalForce('0.5*k*((z-z0)^2)')
        fb_force.addGlobalParameter('k', k)
        fb_force.addGlobalParameter('z0', z_center)

        for atom_index in poly_atoms:
            fb_force.addParticle(atom_index, [])

        omm_sys.addForce(fb_force)

        # Step 3: Add anisotropic barostat to see if higher concentrations are possible for polyatomics
        pressure_vector = Vec3(1.01325, 1.01325, 0.0) * atmosphere
        anisotropic_force=MonteCarloAnisotropicBarostat(pressure_vector,300,True,True,False)

        omm_sys.addForce(anisotropic_force)

        # Save each replicate's results in the dictionary
        omm_objects_mod[f"r{r}"] = {
            'topology': omm_top,
            'positions': omm_pos,
            'system': omm_sys
        }
        
    # Save the full dictionary
    output_picklefile = f"{wdir}/omm_modified_{concentration}.pkl"
    with open(output_picklefile, "wb") as f:
        pickle.dump(omm_objects_mod, f)
    print(f"✅ System dictionary for {concentration}m saved as pickle file")

    return omm_objects_mod
def build_system_HP(concentration, repnum, ion1, ion2, wdir, ff, water, salt_dict, center_atom1, center_atom2, k, delta_z, z_center):
    
    sdf_path1=f'structures/{ion1.lower()}.sdf'
    sdf_path2=f'structures/{ion2.lower()}.sdf'

    if concentration % 1 == 0:
        mi1 = f"{concentration:.1f}"
        mi = str(int(concentration))  # Whole number: strip decimal
    else:
        mi = f"{concentration:.1f}".replace('.', '')  # Float: remove the dot
        mi1 = f"{concentration:.1f}"  # used in salt_dict key lookup

    print(f'Molality of maximum concentration = {mi} mol/kg',
            salt_dict[f'{ion1}{ion2}'][f'Molality {mi1} mol/kg'])
    
    input_picklefile=f"{wdir}/omm_build_{mi}.pkl"

    if os.path.exists(input_picklefile):
        print("File already exists. Skipping code.")
        with open(input_picklefile, "rb") as f:
            omm_build = pickle.load(f)
        return modify_omm_harmonic(mi,repnum,input_picklefile,ion1,ion2,wdir,center_atom1,center_atom2,k,delta_z,z_center)
   
    else:
        print("File not found. Running code...")

        omm_builds = {}

        for r in tqdm(range(repnum), desc=f"Building {mi}m systems", leave=False):
            pdb_path = f'structures/{ion1.lower()}{ion2.lower()}_{mi}m_r{r}.pdb'
            POL1 = Molecule.from_file(sdf_path1)
            POL2 = Molecule.from_file(sdf_path2)
            
            off_top = Topology.from_pdb(pdb_path, unique_molecules=[POL1,POL2])

            if water == 'TIP3P':
                inc = ff.create_interchange(
                    topology=off_top,
                    toolkit_registry=GTR,
                    charge_from_molecules=[water_TIP3P]
                )
            else:
                inc = ff.create_interchange(
                    topology=off_top,
                    toolkit_registry=GTR
                )

            inc.box = boxvectors.get_topology_bbox(off_top)

            omm_top = inc.to_openmm_topology()
            omm_pos = openff_to_openmm(inc.get_positions(include_virtual_sites=True))
            omm_sys = inc.to_openmm_system(combine_nonbonded_forces=False, add_constrained_forces=True)

            # Save each replicate's results in the dictionary
            omm_builds[f"r{r}"] = {
                'topology': omm_top,
                'positions': omm_pos,
                'system': omm_sys
            }

        # Save the full dictionary
        
        with open(input_picklefile, "wb") as f:
            pickle.dump(omm_builds, f)
        print(f"✅ System dictionary for {mi1}m saved as pickle file")
        
        return modify_omm_harmonic(mi,repnum,input_picklefile,ion1,ion2,wdir,center_atom1,center_atom2,k,delta_z,z_center)
    
def simulation_exists(base_dir, postfix):
    """Check if both equil and prod simulation folders already exist."""
    equil_dir = os.path.join(base_dir, f"equil_sim_{postfix}")
    prod_dir = os.path.join(base_dir, f"prod_sim_{postfix}")
    return os.path.exists(equil_dir) and os.path.exists(prod_dir)


def copy_file(source_path, destination_path):
    """
    Copies a file from the source path to the destination path.
    """
    try:
        shutil.copy2(source_path, destination_path)
        print(f"✅ File copied from {source_path} to {destination_path}")
    except FileNotFoundError:
        print(f"❌ Source file {source_path} not found.")
    except PermissionError:
        print(f"❌ Permission denied for {source_path} or {destination_path}.")
    except Exception as e:
        print(f"❌ Error copying file: {e}")

def convert_dcd_to_xtc(input_file, output_file):
    """
    Converts a DCD trajectory file to XTC format using mdconvert.
    """
    command = ['mdconvert', '-o', output_file, input_file]
    try:
        subprocess.run(command, check=True, capture_output=True, text=True)
        print(f"🔁 Successfully converted {input_file} to {output_file}")
    except subprocess.CalledProcessError as e:
        print(f"❌ Error during conversion: {e}")
        print(e.stderr)

def copy_simulation_outputs(molalities, N_replicates, wdir, rdir):
    """
    Copies simulation output files (.pdb) and converts .dcd to .xtc for each molality and replicate.
    """
    wdir = os.path.abspath(wdir)

    for mol in molalities:
        if mol % 1 == 0:
            mi1 = f"{mol:.1f}"
            mi = str(int(mol))  # Whole number: strip decimal
        else:
            mi = f"{mol:.1f}".replace('.', '')  # Float: remove the dot
            mi1 = f"{mol:.1f}"  # used in salt_dict key lookup
        for i in range(N_replicates):
            postfix = f'{mi}m_r{i}'
            pdb_name = f'md{postfix}.pdb'
            xtc_name = f'md{postfix}.xtc'

            # Paths
            source_pdb = os.path.join(wdir, f'prod_sim_{postfix}', f'prod_sim_{postfix}_topology.pdb')
            input_dcd  = os.path.join(wdir, f'prod_sim_{postfix}', f'prod_sim_{postfix}_trajectory.dcd')
            dest_dir   = os.path.join(os.getcwd(), f'{rdir}_{mi}m')
            dest_pdb   = os.path.join(dest_dir, pdb_name)
            dest_xtc   = os.path.join(dest_dir, xtc_name)

            # Create destination folder
            os.makedirs(dest_dir, exist_ok=True)

            # Copy PDB
            print(f"\n📂 Checking: {source_pdb}")
            if os.path.exists(source_pdb):
                copy_file(source_pdb, dest_pdb)
            else:
                print(f"❌ File '{source_pdb}' does not exist.")

            # Convert trajectory
            print(f"🔍 Looking for DCD file: {input_dcd}")
            if os.path.exists(input_dcd):
                convert_dcd_to_xtc(input_dcd, dest_xtc)
            else:
                print(f"❌ DCD file '{input_dcd}' not found.")

RESTRAINT_TYPES = {
    'FBP' : build_system_FBP,
    'HP' : build_system_HP,
}
RESTRAINT_TYPE_ALIASES = {
    'FBP' : 'flat-bottomed potential',
    'HP' : 'harmonic potential',
}

# read arguments from CLI to pass into OpenMM runner

def parse_args() -> Namespace:
    '''Read user inputs from the command line and preprocess basic parts of them'''
    parser = argparse.ArgumentParser(description='Run OpenMM simulation schedules according to presets of simulation parameters')

    parser.add_argument('-n', '--name', help='Tag to use to refer to working directory and internal files', required=True)
    parser.add_argument('-s', '--salt', help='Comma-separated list of ion names to be used', nargs='+', required=True)
    parser.add_argument('-ci', '--centerion', help='Comma-separated list of ions at center of polytomic to be used', nargs='+', required=True)
    parser.add_argument('-w', '--water', help='Water model to be used', required=True)
    parser.add_argument('-rn', '--repnumber', help='Number of replicates to run', required=True)
    parser.add_argument('-m', '--concentrations', help='List of concentrations to be analyzed', nargs='+', required=True)
    parser.add_argument('-c', '--cwd', help='The current working directory that all created files/directories should be made relative to', type=Path, default=Path.cwd())
    parser.add_argument('-spp', '--sim_param_paths',
        help='One or more serialized simulation presets; \
            \nsimulations built from these presets will be run in the order they are provided here',
        type=Path, nargs='+', required=True
    )
    parser.add_argument('-ff', '--ff_files', help='Force field file or files to be used', nargs='+', required=True)
    parser.add_argument('-r', '--restraint_type', help='The kind of potential restraint to apply to ions in the systems', choices=['HP', 'FBP'], required=True)
    parser.add_argument('-k', help='The spring constant to apply in the harmonic potential restraint') # this should also default to a float, but want to retain the NoneType value to check if this is unset
    parser.add_argument('--z_center', help='', type=float, default=7.2)
    parser.add_argument('--delta_z', help='', type=float, default=2.4)
    parser.add_argument('-du', '--distance_unit', help='The unit convention to use for units of length', choices=['angstrom', 'angstroms', 'nanometer', 'nanometers'], default='nanometer')

    args = parser.parse_args()

    salt_env = os.getenv('SALT')
    # Use environment variable if --salt isn't provided
    if not args.salt or args.salt == [""]:
        if salt_env:
            args.salt = salt_env.split()  # Convert "CS BR" -> ['CS', 'BR']
        else:
            raise ValueError("Salt must be provided via '--salt' argument or 'SALT' environment variable.")
    
    ci_env = os.getenv('CENTERIONS')
    # Use environment variable if --salt isn't provided
    if not args.centerion or args.centerion == [""]:
        if ci_env:
            args.centerion = ci_env.split()  # Convert "CS BR" -> ['CS', 'BR']
        else:
            raise ValueError("Center ions must be provided via '--centerion' argument or 'CENTERIONS' environment variable.")
        
    concs_env = os.getenv('CONCENTRATIONS')
    # Use environment variable if --salt isn't provided
    if not args.concentrations or args.concentrations == [""]:
        if concs_env:
            args.concentrations = concs_env.split()  # Convert "CS BR" -> ['CS', 'BR']
        else:
            raise ValueError("List of concentrations must be provided via '--concentrations' argument or 'CONCENTRATIONS' environment variable.")


    ff_env = os.getenv('FF_FILES')    
    # Use environment variable if --ff_files isn't provided 
    if not args.ff_files or args.ff_files == [""]:
        if ff_env:
            args.ff_files = ff_env.split()
        else:
            raise ValueError("Force field files must be provided via '--ff_files' argument or 'FF_FILES' environment variable.")

    # define directories
    args.working_dir = args.cwd / args.name # Can change this to wherever your files are

    # define Quanitity parameters from restraint
    args.restraint_fn = RESTRAINT_TYPES[args.restraint_type]
    logging.info(f'Using restraint type {args.restraint_type} ({RESTRAINT_TYPE_ALIASES[args.restraint_type]})')

    args.distance_unit = getattr(openmm.unit, args.distance_unit) # convert string name of unit into OpenMM object

    args.delta_z = args.delta_z * args.distance_unit
    args.z_center = args.z_center * args.distance_unit

    # slightly more involved handling of force constants, as this default does differ between restraint implementations
    FORCE_CONST_UNIT = (kilojoule_per_mole/nanometer**2)
    if args.k is not None:
        args.k = float(args.k) * FORCE_CONST_UNIT # enforce Quantity with valid units
    else:
        if args.restraint_type == 'HP':
            args.k = 0.68095403 * FORCE_CONST_UNIT
        elif args.restraint_type == 'FBP':
            args.k = 4184 * FORCE_CONST_UNIT
        
    return args


def main():
    '''Run an OpenMM simulation schedule in a given working directory, with arbitrary PDB starting structure and output directory'''
    '''Ensure that a working directory is set up given a name handle that will be used to refer to it
    Returns the resulting path'''
    args = parse_args()
    ion1=args.salt[0]
    ion2=args.salt[1]
    
    water=args.water

    concentration_list = float_list = [float(x) for x in args.concentrations]
    print(concentration_list)
    repnumber=int(args.repnumber)

    salt_dict = load_salt_info(ion1,ion2)

    wdir = args.working_dir
    wdir.mkdir(exist_ok=True)

    rdir = Path(f'{wdir}/result_files')
    

    
    ff_specifiers = [ 
        FF_DIR_REGISTRY['openforcefields'] / f'{args.ff_files[0]}',
        FF_DIR_REGISTRY['openforcefields'] / f'{args.ff_files[1]}'
    ]

    ff_specifiers.append(Path('/scratch/alpine/bamo6610/osmotic-calculations/custom_forcefields/openff-2.1.0-mg.offxml'))

    # check if all ff files exist
    for ff_path in ff_specifiers:
        if not ff_path.exists():
            raise FileNotFoundError(f"Missing force field file: {ff_path}")

    ff = ForceField(*ff_specifiers)

    omm_modobj_all = {}

    for mol in concentration_list:
        key = f'{mol}m'
        omm_modobj = args.restraint_fn( 
        concentration=mol, 
        repnum=repnumber, 
        ion1=ion1,
        ion2=ion2,
        wdir=wdir,
        ff=ff,
        water=water,
        salt_dict=salt_dict,
        center_atom1=args.centerion[0],
        center_atom2=args.centerion[1],
        k=args.k, 
        delta_z=args.delta_z, 
        z_center=args.z_center
        )
        omm_modobj_all[key] = omm_modobj

        if mol % 1 == 0:
            mi1 = f"{mol:.1f}"
            mi = str(int(mol))  # Whole number: strip decimal
        else:
            mi = f"{mol:.1f}".replace('.', '')  # Float: remove the dot
            mi1 = f"{mol:.1f}"  # used in salt_dict key lookup

        for r in range(repnumber):
            postfix=f'{mi}m_r{r}'

            # Skip if simulations already exist
            if simulation_exists(wdir, postfix):
                print(f"Skipping {postfix}, simulation folders already exist.")
                continue

            schedule : dict[str, SimulationParameters] = {
            f'equil_sim_{postfix}' : SimulationParameters(
                integ_params=IntegratorParameters(
                    time_step=2*femtoseconds,
                    total_time=2*nanosecond,
                    num_samples=100,
                ),
                thermo_params=ThermoParameters(
                    ensemble='NVT',
                    temperature=300 * kelvin,
                ),
                reporter_params=ReporterParameters(
                    traj_ext='dcd',
                ),
            ),
            f'prod_sim_{postfix}' : SimulationParameters(
                integ_params=IntegratorParameters(
                    time_step=2*femtoseconds,
                    total_time=20*nanosecond,
                    num_samples=1000,
                ),
                thermo_params=ThermoParameters(
                    ensemble='NVT',
                    temperature=300 * kelvin,
                ),
                reporter_params=ReporterParameters(),
            ),
        }

            omm_top = omm_modobj_all[f'{mi1}m'][f'r{r}']['topology']
            omm_sys = omm_modobj_all[f'{mi1}m'][f'r{r}']['system']
            omm_pos = omm_modobj_all[f'{mi1}m'][f'r{r}']['positions']

            # apply appropriate restraint
            history = run_simulation_schedule(
                working_dir=wdir,
                schedule=schedule,
                init_top=omm_top,
                init_sys=omm_sys,
                init_pos=omm_pos,
                return_history=True
                )

    copy_simulation_outputs(molalities=concentration_list, N_replicates=repnumber, wdir=wdir, rdir=rdir)

if __name__ == '__main__':
    main()