#!/usr/bin/env bash

timedatectl set-timezone America/New_York
timedatectl status
apt-get update
apt-get install -y git htop
echo "ubuntu:2018neuroscience!" | chpasswd

# wait for vagrant to rsync /vagrant
sleep 60

# unzip example directory structures
tar xfz /vagrant/data-shell.tar.gz -C /vagrant

# create additional users, with usernames and passwords provided in text files
while read iuser ipasswd; do

    # Just print this for debugging.
    printf "\tCreating user: %s with password: %s\n" ${iuser} ${ipasswd}

    # Create the user with adduser and copy contents of the datadir into the user's home
    useradd -m -k datadir -s /bin/bash ${iuser}

    # Assign the password to the user.
    # Password is passed via stdin, *twice* (for confirmation).
    passwd ${iuser} <<< "$ipasswd"$'\n'"$ipasswd"

    # give users ssh access (with the same keys as ubuntu)
    cp -R /home/ubuntu/.ssh/ /home/${iuser}/
    chmod 0700 /home/${iuser}/.ssh/
    chmod 0600 /home/${iuser}/.ssh/authorized_keys
    chown -R ${iuser}:${iuser} /home/${iuser}/.ssh/

    # create .Rprofile for each user
    echo ".libPaths( c( .libPaths(), '/home/ubuntu/R-libs') )" >> /home/${iuser}/.Rprofile
    
    # copy example directory structure
    cp -R /vagrant/data-shell/* /home/${iuser}/
    
    # create treasure hunt directories for each user
    mkdir -p /home/${iuser}/treasure
    cp /vagrant/treasureHunt_v2.pl /home/${iuser}/treasure/
    pushd /home/${iuser}/treasure/
    perl treasureHunt_v2.pl
    popd
    rm /home/${iuser}/treasure/treasureHunt_v2.pl

    # set ownership for all files
    chown -R ${iuser} /home/${iuser}/

done < <(paste /vagrant/logins.txt /vagrant/passwords.txt)

# cleanup
rm -rf /vagrant/*

# log file locations
# /var/log/cloud-init.log
# /var/log/cloud-init-output.log

