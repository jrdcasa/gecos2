import os
import utils
import gecos
from pathlib import Path

v_filename = '/home/jramos/Programacion/GITHUB_REPO_DIR/gecos2/examples/09-FolderWithPDBs/PDBs/Conf_Cluster_0001.pdb'
v_nameserver = 'aoki.iem.csic.es'
v_username = 'jramos'
v_keysshfile = '/home/jramos/.ssh/id_rsa_aoki'
v_encrypt_pass = None
v_slurm_part = 'all'
v_bash_extrainfo = ['g16legacy_root=/optnfs/gaussian/g16_legacy','GAUSS_SCRDIR=$TMPDIR','source $g16legacy_root/bsd/g16.profile','export g16legacy_root GAUSS_SCRDIR',]
v_list_nodes = ['']
v_slurm_part_master = 'all'
v_node_master = None
v_maxjobsslurm = 10
v_timelimitslurm = '150:00:00'
v_localdir = '/home/jramos/Programacion/GITHUB_REPO_DIR/gecos2/examples/09-FolderWithPDBs'
v_remotedir = '/optnfs/home/admin_old/jramos/GECOS/99-TEST_PDBFOLDER'
v_pattern = 'Test'
v_databasefullpath = '/home/jramos/Programacion/GITHUB_REPO_DIR/gecos2/examples/09-FolderWithPDBs/Test.db'
v_fileoutputfullpath = '/home/jramos/Programacion/GITHUB_REPO_DIR/gecos2/examples/09-FolderWithPDBs/Test.log'
v_g16path = '/optnfs/gaussian/g16_legacy/g16'
v_g16_keywords = '#p m062x-6-31G* opt'
v_g16_extrainfo = ['',]
v_g16_ncpus = 4
v_g16_mem = 2000
v_g16_charge = 0
v_g16_multiplicity = 1
v_write_gaussian = True
v_run_gaussian = True
v_cutoff_rmsd_qm = 0.0
v_energy_threshold = 9999.9
v_bond_perception = False
v_dockrmsdpack = '/home/jramos/Programacion/sandboxes/sandbox_gecos2/lib/python3.10/site-packages/gecos-2.0-py3.10-linux-x86_64.egg/thirdparty/dockrmsd.x'
v_confpack = 'folder with PDB files'
v_pdbfolder = '/home/jramos/Programacion/GITHUB_REPO_DIR/gecos2/examples/09-FolderWithPDBs/PDBs'

cwd = Path(__file__).parent.absolute()
donepath = os.path.join(cwd, 'done')

if not os.path.isfile(v_databasefullpath) and  not os.path.isfile(donepath):

    log = utils.init_logger(
        "Output2",
        fileoutput=v_fileoutputfullpath,
        append=False, inscreen=True)

    g1 = gecos.GecosPdbFolder(
                           v_pdbfolder,
                           logger=log)

    g1.generate_pdbfolder(
            write_gaussian=v_write_gaussian,
            g16_key=v_g16_keywords,
            g16_nproc=v_g16_ncpus,
            g16_mem=v_g16_mem,
            g16_extra_info=v_g16_extrainfo,
            pattern=v_pattern,
            charge=v_g16_charge,
            multiplicity=v_g16_multiplicity)

    if v_run_gaussian:
        gecos.send_qm_conformers(
            v_nameserver,
            v_databasefullpath,
            v_username,
            v_keysshfile,
            v_localdir,
            v_remotedir,
            v_g16path,
            maxjobsslurm=v_maxjobsslurm,
            regex='*g16*/*.com',
            partition=v_slurm_part,
            exclude_nodes=v_list_nodes,
            ncpus=v_g16_ncpus, 
            partitionmaster=v_slurm_part_master,
            nodemaster=v_node_master,
            mem=v_g16_mem,
            timelimit=v_timelimitslurm,
            encrypted_pass=v_encrypt_pass,
            extraslurminfo=v_bash_extrainfo,
            logger=log)
    else:
        m = 'QM calculations will not be performed'
        print(m) if log is None else log.info(m)
        with open('done', 'w') as f:
            f.writelines('')

else:

    log = utils.init_logger(
            "Output2",
            fileoutput=v_fileoutputfullpath,
            append=True,
            inscreen=True)

    v_outdir = os.path.join(v_localdir, v_pattern + '_g16_results')

    if not os.path.isfile(donepath):
        gecos.check_qm_jobs(
            v_nameserver,
            v_databasefullpath,
            v_username,
            v_keysshfile,
            v_localdir,
            v_remotedir,
            v_outdir,
            v_pattern,
            v_dockrmsdpack,
            encrypted_pass=v_encrypt_pass,
            cutoff_rmsd=v_cutoff_rmsd_qm,
            energy_threshold=v_energy_threshold,
            logger=log)

print("Job Done!!!")
