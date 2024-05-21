import os
import utils
import gecos

v_filename = '/home/jramos/Programacion/GITHUB_REPO_DIR/gecos2/data/IsoP.pdb'
v_nameserver = 'localhost'
v_username = 'jramos'
v_keysshfile = '/home/jramos/.ssh/id_rsa_localhost'
v_encrypt_pass = None
v_slurm_part = 'cpu'
v_list_nodes = ['']
v_slurm_part_master = 'cpu'
v_node_master = None
v_localdir = '/home/jramos/Programacion/GITHUB_REPO_DIR/gecos2/examples/03-IsoP_localhost_rdkit_pm3'
v_remotedir = '/home/jramos/Programacion/GITHUB_REPO_DIR/gecos2/examples/03-IsoP_localhost_rdkit_pm3/remote_files'
v_pattern = '03-IsoP'
v_databasefullpath = '/home/jramos/Programacion/GITHUB_REPO_DIR/gecos2/examples/03-IsoP_localhost_rdkit_pm3/03-IsoP_localhost.db'
v_fileoutputfullpath = '/home/jramos/Programacion/GITHUB_REPO_DIR/gecos2/examples/03-IsoP_localhost_rdkit_pm3/03-IsoP.log'
v_g16path = '/opt/g16/g16'
v_g16_keywords = '#p opt pm3'
v_ncpus = 1
v_mem = 1000
v_charge = 0
v_multiplicity = 1
v_write_gaussian = True
v_nconfs = 10
v_min_iter_mm = 3000
v_cutoff_rmsd_qm = 1.0
v_bond_perception = True
v_dockrmsdpack = '/home/jramos/Programacion/sandboxes/sandbox_gecos2/lib/python3.8/site-packages/gecos-2.0-py3.8-linux-x86_64.egg/thirdparty/dockrmsd.x'
v_confpack = 'rdkit'
v_rdkit_maxattempts = 1000
v_rdkit_prunermsthresh = -0.010
v_rdkit_useexptorsionangleprefs = True
v_rdkit_usebasicknowlwdge = True
v_rdkit_enforcechirality = True
v_rdkit_cluster_method = 'RMSD'
v_rdkit_ffname = 'MMFF'
v_rdkit_cluster_thres = 1.0

if not os.path.isfile(v_databasefullpath) and  not os.path.isfile('done'):

    log = utils.init_logger(
        "Output2",
        fileoutput=v_fileoutputfullpath,
        append=False, inscreen=True)

    g1 = gecos.GecosRdkit(
        filename=v_filename,
        total_charge=v_charge,
        bond_perception=v_bond_perception,
        logger=log)

    g1.generate_conformers(
        v_localdir,
        nconfs=v_nconfs,
        minimize_iterations=v_min_iter_mm,
        maxattempts=v_rdkit_maxattempts,
        prunermsthresh=v_rdkit_prunermsthresh,
        useexptorsionangleprefs=v_rdkit_useexptorsionangleprefs,
        usebasicknowledge=v_rdkit_usebasicknowlwdge,
        enforcechirality=v_rdkit_enforcechirality,
        ff_name=v_rdkit_ffname,
        cluster_method=v_rdkit_cluster_method,
        cluster_threshold=v_rdkit_cluster_thres,
        write_gaussian=v_write_gaussian,
        pattern=v_pattern,
        g16_key=v_g16_keywords,
        g16_nproc=v_ncpus,
        g16_mem=v_mem,
        charge=v_charge,
        multiplicity=v_multiplicity)

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
