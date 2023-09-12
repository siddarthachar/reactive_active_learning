import os
from ase.calculators.orca import ORCA
from ase.calculators.emt import EMT
# MyOrcaProfile = OrcaProfile(["/ihome/crc/install/orca/3.0.3/orca"])
# calc = ORCA(profile=MyOrcaProfile)
from ase import Atoms
from ase.io import read, write
# import ase.calculators.orca as orc

try:
    from ase import Atoms
    import ase.io
    from ase.io import read, write
except ModuleNotFoundError:
    pass

from ase.calculators.lj import LennardJones
from ase.calculators.emt import EMT
from ase.calculators.amber import Amber
from pyGSM.coordinate_systems.delocalized_coordinates import DelocalizedInternalCoordinates
from pyGSM.coordinate_systems.primitive_internals import PrimitiveInternalCoordinates
from pyGSM.coordinate_systems.topology import Topology
from pyGSM.growing_string_methods import DE_GSM
from pyGSM.growing_string_methods import SE_GSM

from pyGSM.level_of_theories.ase import ASELoT
from pyGSM.optimizers.eigenvector_follow import eigenvector_follow
from pyGSM.optimizers.lbfgs import lbfgs
from pyGSM.potential_energy_surfaces import PES
from pyGSM.utilities import nifty
from pyGSM.utilities.elements import ElementData
from pyGSM.molecule import Molecule
from pyGSM.coordinate_systems.slots import Distance, Angle, Dihedral, LinearAngle, MultiAngle
#from deepmd.calculator import DP
# from .main import post_processing, cleanup_scratch

## external isomer init
def get_driving_coord_prim(dc):
    prim = None
    if "ADD" in dc or "BREAK" in dc:
        if dc[1] < dc[2]:
            prim = Distance(dc[1] - 1, dc[2] - 1)
        else:
            prim = Distance(dc[2] - 1, dc[1] - 1)
    elif "ANGLE" in dc:
        if dc[1] < dc[3]:
            prim = Angle(dc[1] - 1, dc[2] - 1, dc[3] - 1)
        else:
            prim = Angle(dc[3] - 1, dc[2] - 1, dc[1] - 1)
    elif "TORSION" in dc:
        if dc[1] < dc[4]:
            prim = Dihedral(dc[1] - 1, dc[2] - 1, dc[3] - 1, dc[4] - 1)
        else:
            prim = Dihedral(dc[4] - 1, dc[3] - 1, dc[2] - 1, dc[1] - 1)
    elif "OOP" in dc:
        # if dc[1]<dc[4]:
        prim = OutOfPlane(dc[1] - 1, dc[2] - 1, dc[3] - 1, dc[4] - 1)
        # else:
        #    prim = OutOfPlane(dc[4]-1,dc[3]-1,dc[2]-1,dc[1]-1)
    return prim

def minimal_wrapper_de_gsm(
        atoms_reactant: Atoms,
#         atoms_product: Atoms, # To make it work for SE-GSM - not needed. 
        calculator,
        fixed_reactant=False,
        fixed_product=False,
):
    # LOT
    # 'EST_Package'
    # 'nproc': args.nproc,
    # 'states': None,

    # PES
    # pes_type = "PES"
    # 'PES_type': args.pes_type,
    # 'adiabatic_index': args.adiabatic_index,
    # 'multiplicity': args.multiplicity,
    # 'FORCE_FILE': args.FORCE_FILE,
    # 'RESTRAINT_FILE': args.RESTRAINT_FILE,
    # 'FORCE': None,
    # 'RESTRAINTS': None,

    # optimizer
    # optimizer_method = "eigenvector_follow"  # OR "lbfgs"
#    optimizer_method = "lbfgs"  # OR "lbfgs"

    # line_search = 'NoLineSearch'  # OR: 'backtrack'
    # only_climb = False
    # 'opt_print_level': args.opt_print_level,
    # step_size_cap = 0.1  # DMAX in the other wrapper

    # molecule
    # coordinate_type = "TRIC" # Allowed: 'Cartesian', 'DLC', 'HDLC', 'TRIC' - what is TRIC? 
    # 'hybrid_coord_idx_file': args.hybrid_coord_idx_file,
    # 'frozen_coord_idx_file': args.frozen_coord_idx_file,
    # 'prim_idx_file': args.prim_idx_file,

    # GSM
    # gsm_type = "SE_GSM"  # SE_GSM, SE_Cross #pbs13
    driving_coordinates=[("ADD", 1, 11), ("ADD", 4, 12)]
    num_nodes = 15  # 20 for SE-GSM #pbs13
    # isomers_file= "isomers.txt",   # driving coordinates, this is a file #pbs13
#     add_node_tol = 0.1  # convergence for adding new nodes
#     conv_tol = 0.0005  # Convergence tolerance for optimizing nodes
    # conv_Ediff = 0.1  # Energy difference convergence of optimization.
    # 'conv_dE': args.conv_dE,
    # conv_gmax = 0.001  # Max grad rms threshold
    # 'BDIST_RATIO': args.BDIST_RATIO,
    # 'DQMAG_MAX': args.DQMAG_MAX,
    # 'growth_direction': args.growth_direction,
    ID = 0
    # 'gsm_print_level': args.gsm_print_level,
    max_gsm_iterations = 100
    max_opt_steps = 30  # 20 for SE-GSM
    # 'use_multiprocessing': args.use_multiprocessing,
    # 'sigma': args.sigma,

    nifty.printcool('Parsed GSM')

    # LOT
    lot = ASELoT.from_options(calculator,
                              geom=[[x.symbol, *x.position] for x in atoms_reactant],
                              ID=ID)

    # PES
    pes_obj = PES.from_options(lot=lot, ad_idx=0, multiplicity=1)

    # Build the topology
    nifty.printcool("Building the topologies")
    element_table = ElementData()
    elements = [element_table.from_symbol(sym) for sym in atoms_reactant.get_chemical_symbols()]

    topology_reactant = Topology.build_topology(
        xyz=atoms_reactant.get_positions(),
        atoms=elements
    )
    
    driving_coord_prims = []
    for dc in driving_coordinates:
        prim = get_driving_coord_prim(dc)
        if prim is not None:
            driving_coord_prims.append(prim)

    for prim in driving_coord_prims:
        if type(prim) == Distance:
            bond = (prim.atoms[0], prim.atoms[1])
            if bond in topology_reactant.edges:
                pass
            elif (bond[1], bond[0]) in topology_reactant.edges():
                pass
            else:
                print(" Adding bond {} to top1".format(bond))
                topology_reactant.add_edge(bond[0], bond[1])
    
#     topology_product = Topology.build_topology(
#         xyz=atoms_product.get_positions(),
#         atoms=elements
#     )

    # Union of bonds
    # debated if needed here or not
#     for bond in topology_product.edges():
#         if bond in topology_reactant.edges() or (bond[1], bond[0]) in topology_reactant.edges():
#             continue
#         print(" Adding bond {} to reactant topology".format(bond))
#         if bond[0] > bond[1]:
#             topology_reactant.add_edge(bond[0], bond[1])
#         else:
#             topology_reactant.add_edge(bond[1], bond[0])

    # primitive internal coordinates
    nifty.printcool("Building Primitive Internal Coordinates")

    addtr = True
    connect = addcart = False
    prim_reactant = PrimitiveInternalCoordinates.from_options(
        xyz=atoms_reactant.get_positions(),
        atoms=elements,
        topology=topology_reactant,
        connect=connect, #coordinate_type == "DLC",
        addtr=addtr, #coordinate_type == "TRIC",
        addcart=addcart, #coordinate_type == "HDLC",
    )

#     prim_product = PrimitiveInternalCoordinates.from_options(
#         xyz=atoms_product.get_positions(),
#         atoms=elements,
#         topology=topology_product,
#         connect=coordinate_type == "DLC",
#         addtr=coordinate_type == "TRIC",
#         addcart=coordinate_type == "HDLC",
#     )

    # add product coords to reactant coords
#     prim_reactant.add_union_primitives(prim_product)

    # Delocalised internal coordinates
    nifty.printcool("Building Delocalized Internal Coordinates")
    deloc_coords_reactant = DelocalizedInternalCoordinates.from_options(
        xyz=atoms_reactant.get_positions(),
        atoms=elements,
        connect=connect, #coordinate_type == "DLC",
        addtr=addtr, #coordinate_type == "TRIC",
        addcart=addcart, #coordinate_type == "HDLC",
        primitives=prim_reactant
    )

    # Molecules
    nifty.printcool("Building the reactant object with TRIC")
    # from_hessian = optimizer_method == "eigenvector_follow"

    molecule_reactant = Molecule.from_options(
        geom=[[x.symbol, *x.position] for x in atoms_reactant],
        PES=pes_obj,
        coord_obj=deloc_coords_reactant,
        Form_Hessian=True,
    )

#     molecule_product = Molecule.copy_from_options(
#         molecule_reactant,
#         xyz=atoms_product.get_positions(),
#         new_node_id=num_nodes - 1,
#         copy_wavefunction=False
#     )

    # optimizer
    nifty.printcool("Building the Optimizer object")
    # opt_options = dict(print_level=1,
    #                    Linesearch=line_search,
    #                    update_hess_in_bg=not (only_climb or optimizer_method == "lbfgs"),
    #                    conv_Ediff=conv_Ediff,
    #                    conv_gmax=conv_gmax,
    #                    DMAX=step_size_cap,
    #                    opt_climb=only_climb)
    # if optimizer_method == "eigenvector_follow":
    #     optimizer_object = eigenvector_follow.from_options(**opt_options)
    # elif optimizer_method == "lbfgs":
    #     optimizer_object = lbfgs.from_options(**opt_options)
    # else:
    #     raise NotImplementedError

    optimizer_object = eigenvector_follow.from_options(Linesearch='backtrack', OPTTHRESH=0.0005, DMAX=0.5, abs_max_step=0.5,
                                                conv_Ediff=0.1)

    # GSM
    nifty.printcool("Building the GSM object")
    gsm = SE_GSM.from_options(
        reactant=molecule_reactant,
#         product=molecule_product,
        nnodes=num_nodes,
#         driving_coords=[("ADD", 11, 13), ("ADD", 11, 14)], #pbs13
        driving_coords=driving_coordinates, #pbs13

#         isomers_file=[("ADD", 4, 12), ("ADD", 1, 11)], #pbs13

#         CONV_TOL=conv_tol,
#         CONV_gmax=conv_gmax,
#         CONV_Ediff=conv_Ediff,
#         ADD_NODE_TOL=add_node_tol,
#         growth_direction=0,  # I am not sure how this works
        optimizer=optimizer_object,
        DQMAG_MAX=0.5, #default value is 0.8
        ADD_NODE_TOL=0.01, #default value is 0.1
        CONV_TOL = 0.0005,
        
#         ID=ID,
#         print_level=1,
#         mp_cores=1,  # parallelism not tested yet with the ASE calculators
#         interp_method="DLC",
    )
    
#     gsm.isomer_init()

# gsm_type = "SE_GSM"  # SE_GSM, SE_Cross #pbs13
#     num_nodes = 20  # 20 for SE-GSM #pbs13
#     isomers_file= "isomers.txt",   # driving coordinates, this is a file #pbs13
#     add_node_tol = 0.1  # convergence for adding new nodes
#     conv_tol = 0.0005  # Convergence tolerance for optimizing nodes
#     conv_Ediff = 100.  # Energy difference convergence of optimization.
#     # 'conv_dE': args.conv_dE,
#     conv_gmax = 100.  # Max grad rms threshold
#     # 'BDIST_RATIO': args.BDIST_RATIO,
#     # 'DQMAG_MAX': args.DQMAG_MAX,
#     # 'growth_direction': args.growth_direction,
#     ID = 0
#     # 'gsm_print_level': args.gsm_print_level,
#     max_gsm_iterations = 100
#     max_opt_steps = 3  # 20 f
#         gsm = SE_GSM.from_options(reactant=M1, nnodes=20, 
#                                   driving_coords=[("ADD", 4, 12), ("ADD", 1, 11)], 
#                                   optimizer=optimizer, print_level=1) #pbs13
    # optimize reactant and product if needed
    
    # SKA31
#     reactant_refE = da_1.get_potential_energy()
    
    if not fixed_reactant:
        nifty.printcool("REACTANT GEOMETRY NOT FIXED!!! OPTIMIZING")
        path = os.path.join(os.getcwd(), 'scratch', f"{ID:03}", "0")
        optimizer_object.optimize(
            molecule=molecule_reactant,
            refE=molecule_reactant.energy,
            opt_steps=100,
            path=path
        )
#     if not fixed_product:
#         nifty.printcool("PRODUCT GEOMETRY NOT FIXED!!! OPTIMIZING")
#         path = os.path.join(os.getcwd(), 'scratch', f"{ID:03}", str(num_nodes - 1))
#         optimizer_object.optimize(
# #             molecule=molecule_product,
#             refE=molecule_product.energy,
#             opt_steps=100,
#             path=path
#         )

    # set 'rtype' as in main one (???)
    # if only_climb:
    #     rtype = 1
    # # elif no_climb:
    # #     rtype = 0
    # else:
    #     rtype = 2

    # do GSM
    nifty.printcool("Main GSM Calculation")
    gsm.go_gsm(max_iters=max_gsm_iterations, opt_steps=max_opt_steps, rtype=2)

    # write the results into an extended xyz file
    string_ase, ts_ase = gsm_to_ase_atoms(gsm)
    ase.io.write(f"opt_converged_{gsm.ID:03d}_ase.xyz", string_ase)
    ase.io.write(f'TSnode_{gsm.ID}.xyz', string_ase)

#     # post processing taken from the main wrapper, plots as well
#     post_processing(gsm, have_TS=True)

#     # cleanup
#     cleanup_scratch(gsm.ID)


def gsm_to_ase_atoms(gsm: SE_GSM):
    # string
    frames = []
    for energy, geom in zip(gsm.energies, gsm.geometries):
        at = Atoms(symbols=[x[0] for x in geom], positions=[x[1:4] for x in geom])
        at.info["energy"] = energy
        frames.append(at)

    # TS
    ts_geom = gsm.nodes[gsm.TSnode].geometry
    ts_atoms = Atoms(symbols=[x[0] for x in ts_geom], positions=[x[1:4] for x in ts_geom])

    return frames, ts_atoms


filepath1 = "/ihome/kjohnson/pbs13/bin/pyGSM/pyGSM/data/diels_alder.xyz"
da_1=read(filepath1, index=0)
# da_2=read(filepath1, index=1)

# minimal_wrapper_de_gsm(da_1, DP(model='/bgfs/kjohnson/ska31/1DeePMD/paper-graphanol/oneside/GOH_DP.pb'),fixed_reactant=True)
minimal_wrapper_de_gsm(da_1, ORCA(charge=0, mult=1, orcasimpleinput='B3LYP 6-31G*',orcablocks='%pal nprocs 1 end\n%scf maxiter 200 end'), fixed_reactant=True)
