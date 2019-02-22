1. Installation

Follow the 'Solr Quick Start' link on the page:

http://lucene.apache.org/solr/resources.html#solr-version-control

to Download/Unpack/Launch Solr on your Solr server machine, if locally, it will be running at:

http://localhost:8983/solr/#/

2. create the core directories
-------------------------------------------------------------------------
1) create folders
mkdir my_solr_install_dir/server/solr/compounds
mkdir -p my_solr_install_dir/server/solr/compounds/conf
mkdir -p my_solr_install_dir/server/solr/compounds/data
mkdir my_solr_install_dir/server/solr/reactions
mkdir -p my_solr_install_dir/server/solr/reactions/conf
mkdir -p my_solr_install_dir/server/solr/reactions/data

-------------------------------------------------------------------------

2) Inside each of the above conf folders, create/copy-paste the solrconfig.xml, schema.xml 
and schema_types.xml.
At this point the data folders are empty.

3) In the browser where you can see the Solr Admin UI, (i.e., http://localhost:8983/solr/#/), 
'Core Admin'-->'Add Core'. Set the instanceDir and dataDir by entering the absolute paths of the dirs created in Step 1):
e.g., for compounds:
        instanceDir= my_solr_install_dir/server/solr/compounds
        dataDir= my_solr_install_dir/server/solr/compounds/data
      for reactions:
        instanceDir= my_solr_install_dir/server/solr/reactions
        dataDir= my_solr_install_dir/server/solr/reactions/data

3. Load data into Solr
At the shell go to where the bin/solr command is and run the following commands to load data
into the corresponding core we have just created in the browser UI. (Note: if you skipped step 3), you will get "HTTP ERROR 404" (Not Found) error.

localhost:my_solr_install_dir qzhang$ bin/post -c core_name path-to-data/datafile.json

When data loading returns error (it happens!), fix the causes and reload until it succeeds.

4. To access the local Solr instance from local modelseed-ui
Inside the config.js of modelseed-ui, add an entry of 'local_solr_url: "http://0.0.0.0:8983/solr/"' to this.services.

Then in file 'modelseed-ui/app/services/biochem.js' and module 'Biochem', find the following line and set the solr_endpoint:

var solr_endpoint = config.services.local_solr_url;

After that, make sure wherever module 'Biochem' is injected and a call to 'Biochem.get*' is changed to'Biochem.get*_solr'.

Then you can see the biochemistry data changes in your local Solr reflected in modelseed-ui.
