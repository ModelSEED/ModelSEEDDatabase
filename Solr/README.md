-------------------------------------------------------------------------
1. Installation of Solr

(I'm assuming that you're a sysadmin, I've not tested how you can install/run solr without sudo access)

Follow the 'Solr Quick Start' link on the page:

http://lucene.apache.org/solr/guide/8_4/solr-tutorial.html

to Download/Unpack/Launch Solr on your Solr server machine, if locally, it will be running at:

http://localhost:8983/solr/#/

As of 02/13/20, I got Solr 8.4.1 running using JDK 11 on Ubuntu 18.04
and found this guide helpful:

https://tecadmin.net/install-apache-solr-on-ubuntu/

As of 04/30/20, I was able to get Solr 8.5 working on my Macbook using MacPorts

On ubuntu, $SOLR_PATH was set to /opt/solr/ and on os x, using MacPorts, $SOLR_PATH was set to /opt/local/share/java/solr-8.5.0/

-------------------------------------------------------------------------
2. Create Solr cores

sudo -u solr ${SOLR_PATH}/bin/solr create_core -c compounds -d cores/compounds/
sudo -u solr ${SOLR_PATH}/bin/solr create_core -c reactions -d cores/reactions/

-------------------------------------------------------------------------
3. Load data into Solr cores

${SOLR_PATH}/bin/post -c compounds ../Biochemistry/compounds.json
${SOLR_PATH}/bin/post -c compounds ../Biochemistry/reactions.json

4. To access the local Solr instance from locally installed modelseed-ui
-------------------------------------------------------------------------
You have to edit the config.js file in the root of the ModelSEED-UI repository, and 
change the solr_url variable from https://modelseed.org/solr/ to http://localhost:8983/solr

NB: Though I don't use it, I found this guide helpful to test basic authentication for our solr endpoint on nginx:
https://docs.nginx.com/nginx/admin-guide/security-controls/configuring-http-basic-authentication/
