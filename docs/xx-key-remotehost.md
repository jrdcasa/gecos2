Tutorial for creating a passwordless RSA key for an SSH connection on a remote computer.
===================================================

#### Create the RSA key

```bash
jramos@totem2:RSA-KEYS$ ssh-keygen -t rsa
Generating public/private rsa key pair.
Enter file in which to save the key (/home/jramos/.ssh/id_rsa): <A NAME FOR THE FILE>
Enter passphrase (empty for no passphrase): 
Enter same passphrase again: 
Your identification has been saved in id_rsa_drago
Your public key has been saved in id_rsa_drago.pub
The key fingerprint is:
SHA256:QDC1BBAoR64T5mePT163HACMdGrzM0kakddgHrS4m0c jramos@totem2
The key's randomart image is:
+---[RSA 3072]----+
| o+=*X*          |
|o.o.X=oo         |
|oo.B.*o          |
|oo. B o.         |
|o. = E .S        |
| .o * o .        |
|   + + . o       |
|    = . o o      |
|     o   o       |
+----[SHA256]-----+


jramos@totem:.ssh$ cat id_rsa_localhost.pub >> ~/.ssh/authorized_keys 

jramos@totem:.ssh$ chmod og-wx ~/.ssh/authorized_keys 
```

This command produce two files:

    id_rsa_<name>
    id_rsa_<name>.pub

```bash
# Local host
scp id_rsa_<name>.pub <username>@<ip>:~/.ssh
# Connect to remote host
ssh <username>@<ip>
cd ~/.ssh
cat id_rsa_<name>.pub >> authorized_keys 

```

Make sure the permissions are:
```bash
chmod g-w /home/user
chmod 700 /home/user/.ssh
chmod 600 /home/user/.ssh/authorized_keys
```

#### Test connection
```bash
jramos@totem:~$ ssh -i .ssh/id_rsa_<name> <username>@<ip>

```

If succesful you connect to remote host without password. 