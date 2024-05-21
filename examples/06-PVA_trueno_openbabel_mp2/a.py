import os
import utils
import gecos

v_filename = '/home/jramos/PycharmProjects/GeCos/examples/06-PVA_trueno_openbabel_mp2/PVA_C5.pdb'
v_nameserver = 'trueno.csic.es'
v_username = 'jramos'
v_keysshfile = '/home/jramos/.ssh/id_rsa_chiripa'
v_encrypt_pass = None
v_slurm_part = 'simacro'
v_list_nodes = ['trueno36', ' trueno37', ' trueno38', ' trueno59']
v_slurm_part_master = 'simacro'
v_node_master = 'trueno36'
v_localdir = '/home/jramos/PycharmProjects/GeCos/examples/06-PVA_trueno_openbabel_mp2/'
v_remotedir = '/home/cfmac/jramos/GECOS/test-06_PVA_openbabel'
v_pattern = '06-PVA'
v_databasefullpath = '/home/jramos/PycharmProjects/GeCos/examples/06-PVA_trueno_openbabel_mp2/06-PVA_trueno.db'
v_fileoutputfullpath = '/home/jramos/PycharmProjects/GeCos/examples/06-PVA_trueno_openbabel_mp2/06-PVA.log'
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
v_confpack = 'openbabel'
v_openbabel_rmsd_cutoff_confab = 0.500000
v_openbabel_energy_cutoff_confab = 50.000000
v_openbabel_verbose = False
v_openbabel_rmsddock_confab = 1.000000
v_openbabel_ffname = 'MMFF'
v_openbabel_cluster_energy_threshold = 99999.000000
v_openbabel_cluster_max_number_cluster = 500

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
