Reference Guide for the Graphical User Interface (GUI) for GeCos
===================================================

#### Molecule input

* **Molecule file (str)**: This field must contain a path to a reference molecule file on your local computer. The pdb file should contain the  section "CONECT". Allowed formats are: *pdb*
* **Name server (str)**: The name or the IP of the remote server. The keyword **localhost** can be used to perform the calculations on local computer. The SLURM system must be installed and configured on either the local or remote servers.
* **Username**: The username for SSH connections to the server.
* **Key SSH file**: Passwordless private RSA key for connecting to the server. (Note 1)
* **Encrypted passwd file**: Encrypted password using the library TODO. (Note 1)
* **SLURM partition**: Name of the partition to which the QM calculation is to be sent.
* **SLURM partition master**: Name of the partition to which the driver (a bash script) that controls the execution of the QM calculations is sent. (Note 2)