-------------------------------------------------------------------------
1. Installation of Solr

Follow the 'Solr Quick Start' link on the page:

http://lucene.apache.org/solr/guide/8_4/solr-tutorial.html

to Download/Unpack/Launch Solr on your Solr server machine, if locally, it will be running at:

http://localhost:8983/solr/#/

NB: As of 02/13/20, we're running Solr 8.4.1 using JDK 11 on Ubuntu 18.04
and we found the first three steps in this guide helpful:

https://tecadmin.net/install-apache-solr-on-ubuntu/

NB: Though we no longer use it, we found this guide helpful to set up basic authentication for our solr endpoint on nginx:
https://docs.nginx.com/nginx/admin-guide/security-controls/configuring-http-basic-authentication/

-------------------------------------------------------------------------
2. Create Solr cores

sudo -u solr /opt/solr/bin/solr create_core -c compounds -d cores/compounds/
sudo -u solr /opt/solr/bin/solr create_core -c reactions -d cores/reactions/

-------------------------------------------------------------------------
3. Load data into Solr cores

/opt/solr/bin/post -u solr -c compounds ../Biochemistry/compounds.json
/opt/solr/bin/post -u solr -c compounds ../Biochemistry/reactions.json

4. To access the local Solr instance from locally installed modelseed-ui
-------------------------------------------------------------------------
Navigate to the root of the ModelSEED-UI repository, and edit the config.js file:

    Change the solr_url from https://modelseed.org/solr/ to http://localhost:8983/solr
