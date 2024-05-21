Tutorial for creating a passwordless RSA key for an SSH connection on a local computer.
===================================================

#### Create the RSA key

```bash
jramos@totem:files$ cd ~/.ssh

jramos@totem:.ssh$ ssh-keygen -t rsa 

Generating public/private rsa key pair.
Enter file in which to save the key (/home/jramos/.ssh/id_rsa): id_rsa_localhost
Enter passphrase (empty for no passphrase): 
Enter same passphrase again: 
Your identification has been saved in id_rsa_localhost
Your public key has been saved in id_rsa_localhost.pub
The key fingerprint is:
SHA256:dR0llTPXo7PSX3WKgbAahHF3DN6F8rtBW6s+Ji0XBgk jramos@totem
The key's randomart image is:
+---[RSA 3072]----+
|    .o. oo... oo=|
|    ..Eoooo. . Bo|
|     . ..*o.. o =|
|      . +.+.oo  o|
|       oSo +.+o.o|
|      .   *.oo. .|
|         o =. . .|
|        o B    . |
|         *..     |
+----[SHA256]-----+

jramos@totem:.ssh$ cat id_rsa_localhost.pub >> ~/.ssh/authorized_keys 

jramos@totem:.ssh$ chmod og-wx ~/.ssh/authorized_keys 
```

Make sure the permissions are:

```bash
chmod g-w /home/user
chmod 700 /home/user/.ssh
chmod 600 /home/user/.ssh/authorized_keys
```

#### Test connection
```bash
jramos@totem:~$ ssh -i .ssh/id_rsa_localhost jramos@totem

 * Documentation:  https://help.ubuntu.com
 * Management:     https://landscape.canonical.com
 * Support:        https://ubuntu.com/advantage
Last login: Sun Mar 20 14:20:33 2022 from 127.0.0.1
jramos@totem:~$ 
```

If succesful you connect to localhost without password. 