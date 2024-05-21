import os
import utils
import gecos
from pathlib import Path

v_filename = '/home/jramos/Programacion/GITHUB_REPO_DIR/gecos2/examples/08-Extract_molecules_fromtrj/frame_101.pdb'
v_nameserver = '161.111.25.96'
v_username = 'jramos'
v_keysshfile = '/home/jramos/.ssh/id_rsa_aoki'
v_encrypt_pass = None
v_slurm_part = 'all'
v_bash_extrainfo = ['hostname','g16root=/optnfs/gaussian','source $g16root/g16_legacy/bsd/g16.profile',]
v_list_nodes = ['']
v_slurm_part_master = 'all'
v_node_master = None
v_maxjobsslurm = 40
v_timelimitslurm = '150:00:00'
v_localdir = '/home/jramos/Programacion/GITHUB_REPO_DIR/gecos2/examples/08-Extract_molecules_fromtrj'
v_remotedir = '/optnfs/home/admin_old/jramos/GECOS/99-TEST'
v_pattern = 'MD_extract'
v_databasefullpath = '/home/jramos/Programacion/GITHUB_REPO_DIR/gecos2/examples/08-Extract_molecules_fromtrj/MD_extract.db'
v_fileoutputfullpath = '/home/jramos/Programacion/GITHUB_REPO_DIR/gecos2/examples/08-Extract_molecules_fromtrj/MD_extract.log'
v_g16path = '/optnfs/gaussian/g16_legacy/g16'
v_g16_keywords = '#p opt b3lyp/6-31G**'
v_g16_extrainfo = ['',]
v_g16_ncpus = 8
v_g16_mem = 8000
v_g16_charge = 0
v_g16_multiplicity = 1
v_write_gaussian = True
v_run_gaussian = True
v_cutoff_rmsd_qm = 0.0
v_energy_threshold = 9999.9
v_bond_perception = True
v_dockrmsdpack = '/home/jramos/Programacion/sandboxes/sandbox_gecos2/lib/python3.8/site-packages/gecos-2.0-py3.8-linux-x86_64.egg/thirdparty/dockrmsd.x'
v_confpack = 'extract from md frame'
v_filemolextractfullpath = "/home/jramos/Programacion/GITHUB_REPO_DIR/gecos2/examples/08-Extract_molecules_fromtrj/frame_101.pdb"
v_extract_method = "Sphere_MINATOM"
v_radius_sphere = 3.0
v_extract_monomers = False
v_extract_pairs = True
v_extract_type = "pair"
v_extract_onlydifferentmols = False

cwd = Path(__file__).parent.absolute()
donepath = os.path.join(cwd, 'done')

if not os.path.isfile(v_databasefullpath) and  not os.path.isfile(donepath):

    log = utils.init_logger(
        "Output2",
        fileoutput=v_fileoutputfullpath,
        append=False, inscreen=True)

    g1 = gecos.GecosExtractNeighbors(
                           v_filemolextractfullpath,
                           write_tcl=True,
                           pattern=v_pattern,
                           radius=v_radius_sphere,
                           logger=log)

    g1.extract_conformers(
            v_localdir,
            calc_com=True,
            method=v_extract_method,
            extract_type=v_extract_type,
            write_gaussian=v_write_gaussian,
            g16_key=v_g16_keywords,
            g16_nproc=v_g16_ncpus,
            g16_mem=v_g16_mem,
            g16_extra_info=v_g16_extrainfo,
            pattern=v_pattern,
            charge=v_g16_charge,
            multiplicity=v_g16_multiplicity,
            onlydifferentmols=v_extract_onlydifferentmols)

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
