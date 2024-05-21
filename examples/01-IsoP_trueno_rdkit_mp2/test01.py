import os
import utils
import gecos

v_filename = '/home/jramos/PycharmProjects/GeCos/data/IsoP.pdb'
v_nameserver = 'trueno.csic.es'
v_username = 'jramos'
v_keysshfile = '/home/jramos/.ssh/id_rsa_chiripa'
v_encrypt_pass = 'None'
v_slurm_part = 'simacro'
v_list_nodes = ['trueno36', 'trueno37', 'trueno38', 'trueno59']
v_slurm_part_master = 'simacro'
v_node_master = 'trueno36'
v_localdir = '/home/jramos/PycharmProjects/GeCos/examples/01-IsoP_trueno_rdkit_mp2'
v_remotedir = '/home/cfmac/jramos/GECOS/test-01_IsoP'
v_pattern = '01-IsoP'
v_databasefullpath = '/home/jramos/PycharmProjects/GeCos/examples/01-IsoP_trueno_rdkit_mp2/01-IsoP_trueno.db'
v_fileoutputfullpath = '/home/jramos/PycharmProjects/GeCos/examples/01-IsoP_trueno_rdkit_mp2/01-IsoP.log'
v_g16path = '/opt/gaussian/g16_legacy/'
v_g16_keywords = '#p opt mp2/6-31g(d,p)'
v_ncpus = 1
v_mem = 4000
v_charge = 0
v_multiplicity = 1
v_write_gaussian = True
v_nconfs = 4000
v_min_iter_mm = 3000
v_cutoff_rmsd_qm = 1.0
v_bond_perception = True
v_dockrmsdpack = '/home/jramos/PycharmProjects/GeCos/thirdparty/dock_rmsd/dockrmsd.x'
v_confpack = 'rdkit'
v_rdkit_maxattempts = 1000
v_rdkit_prunermsthresh = -0.01
v_rdkit_useexptorsionangleprefs = True
v_rdkit_usebasicknowlwdge = True
v_rdkit_enforcechirality = True
v_rdkit_cluster_method = 'rmsd'
v_rdkit_ffname = 'MMFF'
v_rdkit_cluster_thres = 1.0

if not os.path.isfile(v_databasefullpath):

    log = utils.init_logger(
        "Output2",
        fileoutput=v_fileoutputfullpath,
        append=False, inscreen=False)

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

    log = utils.init_logger(
            "Output2",
            fileoutput=v_fileoutputfullpath,
            append=True,
            inscreen=False)

    v_outdir = os.path.join(v_localdir, v_pattern + '_g16_conformers')

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
