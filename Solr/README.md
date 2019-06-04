1. Installation

Follow the 'Solr Quick Start' link on the page:

http://lucene.apache.org/solr/resources.html#solr-version-control

to Download/Unpack/Launch Solr on your Solr server machine, if locally, it will be running at:

http://localhost:8983/solr/#/

2. create the core directories
-------------------------------------------------------------------------
e.g., on my laptop suppose Solr is installed under my_solr_install_dir and
my_solr_install_dir=/Users/qzhang/SOLR/solr-7.7.0

2.1) create one folder for each core (and subfolders of conf and data for that core)

    mkdir my_solr_install_dir/server/solr/compounds
    mkdir -p my_solr_install_dir/server/solr/compounds/conf
    mkdir -p my_solr_install_dir/server/solr/compounds/data
    mkdir my_solr_install_dir/server/solr/reactions
    mkdir -p my_solr_install_dir/server/solr/reactions/conf
    mkdir -p my_solr_install_dir/server/solr/reactions/data

-------------------------------------------------------------------------

2.2) Inside each of the above conf folders, create/copy-paste the solrconfig.xml and elevate.xml from
    my_solr_install_dir/server/solr/configsets/sample_techproducts_configs/conf/.

And to keep data consistency for ModelSEEDDatabase, copy-paste files schema.xml and schema_types.xml
from the corresponding dirs of the model_compound and model_reaction of the patricbrc github repo
(https://github.com/PATRIC3/patric_solr) into folders compounds/conf/ and reactions/conf/ respectively.

Because the patricbrc repo was based on Solr5.3, we need to replace 'solr.TextDocValueField' with
'solr.SortableTextField' in file schema_types.xml.

At this point the data folders are empty.

2.3) In the browser where you can see the Solr Admin UI, (i.e., http://localhost:8983/solr/#/), 
Under 'Core Admin'-->'Add Core'. Set the instanceDir and dataDir by entering the absolute paths of the dirs created in Step 2.1):
e.g., for compounds:
        instanceDir= my_solr_install_dir/server/solr/compounds
        dataDir= my_solr_install_dir/server/solr/compounds/data
      for reactions:
        instanceDir= my_solr_install_dir/server/solr/reactions
        dataDir= my_solr_install_dir/server/solr/reactions/data

3. Load data into Solr
At the shell go to folder `my_solr_install_dir/bin/` and run the following commands to load data
into the corresponding core we have just created in the browser UI. [Note: if you skipped step 2.3),
you will get "HTTP ERROR 404" (Not Found) error.]

    localhost:my_solr_install_dir qzhang$ bin/post -c core_name path-to-data/datafile.json

When data loading returns error (it happens, probably due to file format or acceptable string values, etc.), 
fix the causes and reload until it succeeds.

4. To access the local Solr instance from locally installed modelseed-ui
Inside the config.js of modelseed-ui folder, add an entry of 'local_solr_url: "http://0.0.0.0:8983/solr/"' to this.services.

Then in file 'modelseed-ui/app/services/biochem.js' and module 'Biochem', find the following line and set the solr_endpoint:

var solr_endpoint = config.services.local_solr_url;

After that, make sure wherever module 'Biochem' is injected and a call to 'Biochem.get*' is changed to'Biochem.get*_solr'.

Then you can see the biochemistry data changes in your local Solr reflected in modelseed-ui.

-------------------------------------------------------------------------
5. The above procedures were tested for running Solr on my Mac laptop. They should apply to setting up Solr on the megilen server,
where the cloud version of Solr (7) has been installed.

For example, run the following command to load the data from compounds.json in the `compounds` core

    megilen-cloud-server:modelseedSolrUser$ megilen_solr_installation_dir/bin/post -c compounds path-to-biochem_data/compounds.json
