import os
import utils
import gecos

v_filename = '/home/jramos/Programacion/GITHUB_REPO_DIR/gecos2/examples/08-Extract_molecules/EVOH_residues.pdb'
v_nameserver = 'localhost'
v_username = 'jramos'
v_keysshfile = '/home/jramos/.ssh/id_rsa_localhost'
v_encrypt_pass = None
v_slurm_part = 'cpu'
v_list_nodes = ['']
v_slurm_part_master = 'cpu'
v_node_master = None
v_localdir = '/home/jramos/Programacion/GITHUB_REPO_DIR/gecos2/examples/08-Extract_molecules'
v_remotedir = '/home/jramos/Programacion/GITHUB_REPO_DIR/gecos2/examples/08-Extract_molecules/remote_dir'
v_pattern = 'MD_extract'
v_databasefullpath = '/home/jramos/Programacion/GITHUB_REPO_DIR/gecos2/examples/08-Extract_molecules/MD_extract.db'
v_fileoutputfullpath = '/home/jramos/Programacion/GITHUB_REPO_DIR/gecos2/examples/08-Extract_molecules/MD_extract.log'
v_g16path = '/opt/g16/g16'
v_g16_keywords = '#p opt pm3'
v_ncpus = 1
v_mem = 1000
v_charge = 0
v_multiplicity = 1
v_write_gaussian = True
v_nconfs = 2000
v_min_iter_mm = 1000
v_cutoff_rmsd_qm = 0.5
v_bond_perception = True
v_dockrmsdpack = '/home/jramos/Programacion/sandboxes/sandbox_gecos2/lib/python3.8/site-packages/gecos-2.0-py3.8-linux-x86_64.egg/thirdparty/dockrmsd.x'
v_confpack = 'extract from md frame'
v_filemolextractfullpath = "/home/jramos/Programacion/GITHUB_REPO_DIR/gecos2/examples/08-Extract_molecules/EVOH_residues.pdb"
v_extract_method = "Sphere"
v_radius_sphere =        7.0
v_extract_monomers = True
v_extract_pairs = False
v_extract_trimers = False

if not os.path.isfile(v_databasefullpath) and  not os.path.isfile('done'):

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
            write_gaussian=v_write_gaussian,
            g16_key=v_g16_keywords,
            g16_nproc=v_ncpus,
            g16_mem=v_mem,
            pattern=v_pattern,
            charge=v_charge,
            multiplicity=v_multiplicity
     )

    if v_write_gaussian:
        gecos.send_qm_conformers(
            v_nameserver,
            v_databasefullpath,
            v_username,
            v_keysshfile,
            v_localdir,
            v_remotedir,
            v_g16path,
            regex='*g16*/*.com',
            partition=v_slurm_part,
            exclude_nodes=v_list_nodes,
            ncpus=v_ncpus, 
            partitionmaster=v_slurm_part_master,
            nodemaster=v_node_master,
            mem=v_mem,
            encrypted_pass=v_encrypt_pass,
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

    v_outdir = os.path.join(v_localdir, v_pattern + '_g16_conformers')

    if not os.path.isfile('done'):
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
            logger=log)

print("Job Done!!!")
