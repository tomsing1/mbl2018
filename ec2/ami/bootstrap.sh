#!/usr/bin/env bash

timedatectl set-timezone America/New_York
timedatectl status

# delay to finish system setup
sleep 60

apt-get update
apt-get install -y git python-pip htop apache2 emacs nvme-cli
echo "ubuntu:2018neuroscience!" | chpasswd  # update default password!

# install AWS CLI tool
pip install --upgrade pip
pip install --upgrade awscli

# install conda
wget -q https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh \
  -O miniconda.sh
chmod +x miniconda.sh
./miniconda.sh -b -p /home/ubuntu/miniconda3
chown -R ubuntu /home/ubuntu/
echo PATH=/home/ubuntu/miniconda3/bin:$PATH >> /etc/environment

# configure apache2 webserver to listen on port 8080 and show ubuntu's home directory
cat <<EOF >> /etc/apache2/apache2.conf
 <Directory /home/ubuntu/>
       Options Indexes FollowSymLinks
       AllowOverride None
       Require all granted
</Directory>
EOF
echo "Listen 8080" > /etc/apache2/ports.conf
sed -i s/'DocumentRoot \/var\/www\/html'/'DocumentRoot \/home\/ubuntu'/ \
    /etc/apache2/sites-available/000-default.conf
sed -i s/'<VirtualHost \*:80>'/'<VirtualHost \*:8080>'/ \
    /etc/apache2/sites-available/000-default.conf
/etc/init.d/apache2 restart

# clean AMI
apt-get clean
rm -rf /root/.bash_history /root/.cache /root/.viminfo /root/startup.log \
  /usr/local/Rmpi/lastdeviceused /root/display.txt /root/xvfb.log \
  /root/.subversion /usr/local/etc/hostfile.txt /usr/local/lib/R/etc/Renviron.site \
  /root/.aws
rm -rf /home/ubuntu/.bash_history /home/ubuntu/.sudo_as_admin_successful \
  /home/ubuntu/.cache /home/ubuntu/.viminfo /home/ubuntu/.Rhistory \
  /home/ubuntu/.rstudio /home/ubuntu/.subversion /home/ubuntu/.mpi_is_set_up \
  /home/ubuntu/.aws
history -c
su -c "history -c" -s /bin/bash ubuntu

# log file: /var/log/cloud-init-output.log
