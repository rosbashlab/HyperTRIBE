MariaDB Set Up
==============

Source of `sample code <https://www.vultr.com/docs/how-to-install-apache-24-mariadb-10-and-php-7x-on-ubuntu-16-04>`_

Install MariaDB 10.x
--------------------

You can use the following commands to install MariaDB 10.1 on your Ubuntu 16.04 x64 system.

Setup the system apt repo:
::

    sudo apt-get install software-properties-common
    sudo apt-key adv --recv-keys --keyserver hkp://keyserver.ubuntu.com:80 0xF1656F24C74CD1D8
    sudo add-apt-repository 'deb [arch=amd64,i386,ppc64el] http://mirror.jmu.edu/pub/mariadb/repo/10.1/ubuntu xenial main'

Install MariaDB:
::

    sudo apt update -y
    sudo apt install -y mariadb-server

During the installation process, the MariaDB package configuration wizard will automatically pop up and ask you to setup a new password for the MariaDB **root** user. For now, just press **Enter** every time the wizard pops up to skip this step because we will setup a password for the MariaDB **root** user in the following securing MariaDB procedure.

Having MariaDB installed, you can confirm the installation with:
::

    mysql -V

The output should be similar to:
::

    mysql  Ver 15.1 Distrib 10.1.22-MariaDB, for debian-linux-gnu (x86_64) using readline 5.2

Start the MariaDB service:
::

    sudo systemctl start mariadb.service
    sudo systemctl enable mariadb.service

Secure the installation of MariaDB:
::

    sudo /usr/bin/mysql_secure_installation

During the interactive process, answer questions one by one as follows:
::

    Enter current password for root (enter for none): <Enter>
    Set root password? [Y/n]: Y
    New password: <your-MariaDB-root-password>
    Re-enter new password: <your-MariaDB-root-password>
    Remove anonymous users? [Y/n]: Y
    Disallow root login remotely? [Y/n]: Y
    Remove test database and access to it? [Y/n]: Y
    Reload privilege tables now? [Y/n]: Y

Note: Be sure to replace *<your-MariaDB-root-password>* with your own MariaDB root password.

In this fashion, MariaDB 10.1 has been securely installed onto your system. In the future, you can setup designated users and databases for your web apps as follows:

Log into the MySQL shell as *root*:
::

    mysql -h localhost -u root -p

Type the MariaDB root password you set earlier when prompted.


Create a MariaDB database 
::
    CREATE DATABASE dmseq;

Feel free to create a database with a different name, just ensure that the same  database is used for as well.

