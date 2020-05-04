# ModelSEEDDatabase in Solr

The primary purpose of this directory is to store and maintain the
schema that we use for our ModelSEEDDatabase Solr cores.

The same files can be used in anyone's local instance of Solr,
especially if they wish to maintain their own version of the
database. So I explain here how to do this.

(The assumption is that the user is a sysadmin or has sudo access,
I've not attempted to install/run solr without sudo access)

## Installation of Solr

Follow the 'Solr Quick Start' link on the page:  
http://lucene.apache.org/solr/guide/8_4/solr-tutorial.html

to Download/Unpack/Launch Solr on your Solr server machine, if
locally, it will be running at:  http://localhost:8983/solr/

As of 02/13/20, I got Solr 8.4.1 running using JDK 11 on Ubuntu 18.04
and found this guide helpful:  
https://tecadmin.net/install-apache-solr-on-ubuntu/

As of 04/30/20, I was able to get Solr 8.5 working on my Macbook using
MacPorts:  https://www.macports.org/

For the purpose of this guide, I'm using the `$SOLR_PATH` environment
variable, but it differed between the two environments. On Ubuntu,
`$SOLR_PATH` is `/opt/solr/` and on OS X, using MacPorts, `$SOLR_PATH`
is `/opt/local/share/java/solr-8.5.0/`.

## Create Solr Cores

Here we create the cores using the schema that is available in this directory

```
sudo -u solr ${SOLR_PATH}/bin/solr create_core -c compounds -d cores/compounds/
sudo -u solr ${SOLR_PATH}/bin/solr create_core -c reactions -d cores/reactions/
```

## Load Biochemistry Data

Here we load the actual data:

```
${SOLR_PATH}/bin/post -c compounds ../Biochemistry/compounds.json
${SOLR_PATH}/bin/post -c compounds ../Biochemistry/reactions.json
```

## Access the Data

You can query the data in a number of ways:

1. Solr's default UI:  
http://localhost:8983/solr/#/compounds/query  
http://localhost:8983/solr/#/reactions/query

2. ModelSEED-UI:  
The ModelSEED has a custom search interface for exploring the
biochemistry:  
https://modelseed.org/biochem/reactions/<p>The code
for the web server is publicly available:  
https://github.com/ModelSEED/ModelSEED-UI</p><p>and any user can use
this code to host their own ModelSEED website. In doing so, they can
point it to their local Solr instance, and use it to explore their
biochemical data. To do so, you need to update a line in this file:  
https://github.com/ModelSEED/ModelSEED-UI/blob/master/config.js</p><p>You
have to change the variable `solr_url` so that its configured to point
to  
http://localhost:8983/solr</p>

3. <a href="https://en.wikipedia.org/wiki/Representational_state_transfer">REST</a>  
Apache Solr is RESTful, meaning it has an API that allows the server
to be queried using `http` and any command-line tool or programming
language. We have examples of this, in using `curl` and `python`
in the _examples_ folder:  
```
examples/Check_Status_of_Cores.sh
examples/Check_Status_of_Cores.py
examples/Query_ModelSEED_Example.py
examples/Query_ModelSEED_Example.sh
```

## KBase

The same biochemistry data hosted here and at https://modelseed.org is
also available at KBase but the data is hosted by a different
database. KBase has a python client that can be used to access the
data, and the instructions for doing so are in the example code in the same folder:

```
examples/Query_KBase_Example.py
```

## Authentication

We're using our own approach to restrict access to the Solr
configuration, and so we are not using basic authentication, but I
found this guide helpful for testing solr authentication:

https://docs.nginx.com/nginx/admin-guide/security-controls/configuring-http-basic-authentication/
