1. allow ssh login to localhost

a. generate key pair    ssh-keygen 
b. cat .ssh/id_rsa.pub >> authorized_keys
c. test connection: ssh admin@localhost  --> should log in.  after saying "Yes" to question


2. install snow package

sudo R CMD INSTALL mlds/CRAN/snow_0.3-13.tar.gz


3. run nosetests
