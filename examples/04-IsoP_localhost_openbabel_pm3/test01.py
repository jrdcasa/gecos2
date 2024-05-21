import os
import utils
import gecos

v_filename = '/home/jramos/Programacion/GITHUB_REPO_DIR/gecos2/data/IsoP.pdb'
v_nameserver = 'localhost'
v_username = 'jramos'
v_keysshfile = '/home/jramos/Programacion/GITHUB_REPO_DIR/gecos2/examples/03-IsoP_localhost_rdkit_pm3/keys/keylocal.pem'
v_encrypt_pass = '/home/jramos/Programacion/GITHUB_REPO_DIR/gecos2/examples/03-IsoP_localhost_rdkit_pm3/keys/passwd_encrypted.bin'
v_slurm_part = 'cpu'
v_list_nodes = ['']
v_slurm_part_master = 'cpu'
v_node_master = 'totem'
v_localdir = '/home/jramos/Programacion/GITHUB_REPO_DIR/gecos2/examples/04-IsoP_localhost_openbabel_pm3'
v_remotedir = '/home/jramos/Programacion/GITHUB_REPO_DIR/gecos2/examples/04-IsoP_localhost_openbabel_pm3/remote_files'
v_pattern = '04-IsoP'
v_databasefullpath = '/home/jramos/Programacion/GITHUB_REPO_DIR/gecos2/examples/04-IsoP_localhost_openbabel_pm3/04-IsoP_localhost.db'
v_fileoutputfullpath = '/home/jramos/Programacion/GITHUB_REPO_DIR/gecos2/examples/04-IsoP_localhost_openbabel_pm3/04-IsoP.log'
v_g16path = '/opt/g16/g16'
v_g16_keywords = '#p opt pm3'
v_ncpus = 1
v_mem = 1000
v_charge = 0
v_multiplicity = 1
v_write_gaussian = True
v_nconfs = 500
v_min_iter_mm = 3000
v_cutoff_rmsd_qm = 1.0
v_bond_perception = True
v_dockrmsdpack = '/home/jramos/Programacion/GITHUB_REPO_DIR/gecos2/thirdparty/dock_rmsd/dockrmsd.x'
v_confpack = 'openbabel'
v_openbabel_rmsd_cutoff_confab = 0.5
v_openbabel_energy_cutoff_confab = 50.0
v_openbabel_verbose = False
v_openbabel_rmsddock_confab = 2.0
v_openbabel_ffname = 'MMFF'
v_openbabel_cluster_energy_threshold = 99999.0
v_openbabel_cluster_max_number_cluster = 100

if not os.path.isfile(v_databasefullpath):

    log = utils.init_logger(
        "Output2",
        fileoutput=v_fileoutputfullpath,
        append=False, inscreen=False)

    g1 = gecos.GecosPyBabel(
        filename=v_filename,
        exec_rmsddock=v_dockrmsdpack,
        total_charge=v_charge,
        bond_perception=v_bond_perception,
        logger=log)

    g1.generate_conformers(
        v_localdir,
        nconfs=v_nconfs,
        minimize_iterations=v_min_iter_mm,
        rmsd_cutoff_confab=v_openbabel_rmsd_cutoff_confab,
        energy_cutoff_confab=v_openbabel_energy_cutoff_confab,
        confab_verbose_confab=v_openbabel_verbose,
        cutoff_rmsddock_confab=v_openbabel_rmsddock_confab,
        energy_threshold_cluster=v_openbabel_cluster_energy_threshold,
        max_number_cluster=v_openbabel_cluster_max_number_cluster,
        ff_name=v_openbabel_ffname,
        pattern=v_pattern,
        write_gaussian=v_write_gaussian,
        g16_key=v_g16_keywords,
        g16_nproc=v_ncpus,
        g16_mem=v_mem,
        charge=v_charge,
        multiplicity=v_multiplicity
        )


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
